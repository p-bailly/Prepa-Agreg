// Séance du jour : cartes à revoir puis cartes neuves (plafond quotidien), filtres, ordre d'introduction.

import { cleVolet, rangNiveau, voletsDe, type Carte, type Contenu, type Volet } from '../contenu/types'
import type { Progression } from '../lib/progression'
import type { Reglages } from '../lib/reglages'
import { joursEntre } from '../lib/temps'
import { etatAffiche, etatCarte, prerequisManquants, sousLeDepart, voletsDus } from './etat'

export interface ElementSeance {
  carte: string
  /** Volets à revoir, le premier étant le volet principal (noté par le pied de l'écran). */
  volets: Volet[]
  neuve: boolean
}

export type ReglagesSeance = Pick<
  Reglages,
  'niveauDepart' | 'niveauxInferieurs' | 'niveauxActifs' | 'plafondNeuves' | 'seanceCartes' | 'filtres'
>

/** Durées estimées, pour « environ N minutes ». */
export const MINUTES_REVISION = 0.75
export const MINUTES_NEUVE = 1.5

/**
 * Ordre d'introduction des cartes neuves : niveau le plus bas d'abord ; dans un niveau, les thèmes alternent (première
 * carte de chaque thème, puis deuxième…) et chaque thème suit l'ordre de son fichier.
 */
export function ordreIntroduction(contenu: Contenu): Carte[] {
  const rangTheme = new Map(contenu.themes.map((t, i) => [t.id, i]))
  const position = new Map<string, number>()
  const compteur = new Map<string, number>()
  for (const c of contenu.cartes) {
    const cle = `${c.theme}/${c.niveau}`
    const n = compteur.get(cle) ?? 0
    position.set(c.id, n)
    compteur.set(cle, n + 1)
  }
  return [...contenu.cartes].sort(
    (a, b) =>
      rangNiveau(a.niveau) - rangNiveau(b.niveau) ||
      position.get(a.id)! - position.get(b.id)! ||
      (rangTheme.get(a.theme) ?? 99) - (rangTheme.get(b.theme) ?? 99),
  )
}

/** La carte passe-t-elle les réglages (niveaux actifs, niveau de départ) et les filtres de séance ? */
export function carteRetenue(c: Carte, p: Progression, r: ReglagesSeance): boolean {
  if (!r.niveauxActifs.includes(c.niveau) || sousLeDepart(c, r)) return false
  const f = r.filtres
  if (f.themes.length && !f.themes.includes(c.theme)) return false
  if (f.niveaux.length && !f.niveaux.includes(c.niveau)) return false
  if (f.etats.length && !f.etats.includes(etatAffiche(c, p.volets, r))) return false
  return true
}

/** Cartes neuves introduites aujourd'hui (pour le plafond quotidien). */
export function neuvesDuJour(contenu: Contenu, p: Progression, jour: string): number {
  return contenu.cartes.filter((c) => {
    const premieres = voletsDe(c)
      .map((v) => p.volets[cleVolet(c.id, v)]?.premiere)
      .filter((x): x is string => !!x)
    return premieres.length > 0 && premieres.sort()[0] === jour
  }).length
}

/** Cartes à revoir, les plus en retard d'abord. */
export function cartesARevoir(contenu: Contenu, p: Progression, r: ReglagesSeance, maintenant: number, jour: string): ElementSeance[] {
  const retard = (c: Carte, volets: Volet[]) =>
    Math.max(...volets.map((v) => {
      const e = p.volets[cleVolet(c.id, v)]
      return e ? joursEntre(e.echeance, jour) : 0
    }))
  return contenu.cartes
    .filter((c) => carteRetenue(c, p, r))
    .map((c) => ({ c, volets: voletsDus(c, p.volets, maintenant, jour) }))
    .filter((x) => x.volets.length > 0)
    .map((x) => ({ ...x, retard: retard(x.c, x.volets) }))
    .sort((a, b) => b.retard - a.retard)
    .map((x) => ({ carte: x.c.id, volets: x.volets, neuve: false }))
}

/**
 * Cartes neuves proposables aujourd'hui, dans l'ordre d'introduction, au plus `max`. Une carte neuve n'est proposée que
 * si ses prérequis sont au moins « en cours » ; une carte introduite plus tôt dans la même liste compte comme en cours.
 */
export function cartesNeuves(contenu: Contenu, p: Progression, r: ReglagesSeance, max: number): ElementSeance[] {
  const parId = new Map(contenu.cartes.map((c) => [c.id, c]))
  const choisies: Carte[] = []
  const prises = new Set<string>()
  const candidates = ordreIntroduction(contenu).filter((c) => etatCarte(c, p.volets) === 'neuve' && carteRetenue(c, p, r))
  // Plusieurs passes : une carte débloquée par une carte choisie plus loin dans l'ordre peut entrer à la passe suivante.
  let ajout = true
  while (ajout && choisies.length < max) {
    ajout = false
    for (const c of candidates) {
      if (choisies.length >= max) break
      if (prises.has(c.id) || prerequisManquants(c, parId, p.volets, r, prises).length) continue
      choisies.push(c)
      prises.add(c.id)
      ajout = true
    }
  }
  return choisies.map((c) => ({ carte: c.id, volets: voletsDe(c), neuve: true }))
}

export interface ResumeJour {
  aRevoir: number
  nouvelles: number
  /** Durée estimée de la séance que lance « Réviser ». */
  minutes: number
}

export function plafondRestant(contenu: Contenu, p: Progression, r: ReglagesSeance, jour: string): number {
  return Math.max(0, r.plafondNeuves - neuvesDuJour(contenu, p, jour))
}

export function resumeDuJour(
  contenu: Contenu,
  p: Progression,
  r: ReglagesSeance & Pick<Reglages, 'seanceMinutes'>,
  maintenant: number,
  jour: string,
): ResumeJour {
  const aRevoir = cartesARevoir(contenu, p, r, maintenant, jour).length
  const nouvelles = cartesNeuves(contenu, p, r, plafondRestant(contenu, p, r, jour)).length
  const seance = construireSeance(contenu, p, r, maintenant, jour)
  const estimees = seance.reduce((s, e) => s + (e.neuve ? MINUTES_NEUVE : MINUTES_REVISION), 0)
  return { aRevoir, nouvelles, minutes: Math.min(r.seanceMinutes, Math.ceil(estimees)) }
}

/**
 * Séance courte : au plus `seanceCartes` cartes, les révisions d'abord servies, complétées par des cartes neuves dans la
 * limite du plafond du jour. Les neuves sont réparties entre les révisions, en gardant leur ordre (prérequis d'abord).
 */
export function construireSeance(contenu: Contenu, p: Progression, reglages: ReglagesSeance, maintenant: number, jour: string): ElementSeance[] {
  const revisions = cartesARevoir(contenu, p, reglages, maintenant, jour).slice(0, reglages.seanceCartes)
  const place = Math.min(reglages.seanceCartes - revisions.length, plafondRestant(contenu, p, reglages, jour))
  const neuves = place > 0 ? cartesNeuves(contenu, p, reglages, place) : []
  if (!neuves.length) return revisions
  if (!revisions.length) return neuves
  const total = revisions.length + neuves.length
  // La k-ième neuve se place vers (k + ½) × total / nombre de neuves : positions distinctes, réparties.
  const positions = new Set(neuves.map((_, k) => Math.floor(((k + 0.5) * total) / neuves.length)))
  const sortie: ElementSeance[] = []
  let r = 0
  let n = 0
  for (let i = 0; i < total; i++) {
    if ((positions.has(i) && n < neuves.length) || r >= revisions.length) sortie.push(neuves[n++])
    else sortie.push(revisions[r++])
  }
  return sortie
}
