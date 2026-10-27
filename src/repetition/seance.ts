// Journée : la leçon du jour (notions nouvelles, plafond quotidien) puis la séance de révision (cartes déjà apprises),
// filtres, ordre d'introduction. On apprend d'abord ; une carte n'entre en révision qu'une fois apprise.

import { rangNiveau, type Carte, type Chapitre, type Contenu, type Volet } from '../contenu/types'
import type { Progression } from '../lib/progression'
import type { Reglages } from '../lib/reglages'
import { joursEntre } from '../lib/temps'
import { etatAffiche, etatCarte, introduction, prerequisManquants, sousLeDepart, voletsDus } from './etat'

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

/** Durées estimées, pour « environ N minutes » : une révision, une notion apprise en leçon. */
export const MINUTES_REVISION = 0.75
export const MINUTES_NOTION = 3

/**
 * Ordre d'introduction des cartes neuves : niveau le plus bas d'abord ; dans un niveau, les thèmes alternent (première
 * carte de chaque thème, puis deuxième…) et chaque thème suit l'ordre de ses chapitres (themes.json), puis l'ordre du
 * fichier dans un chapitre.
 */
export function ordreIntroduction(contenu: Contenu): Carte[] {
  const rangTheme = new Map(contenu.themes.map((t, i) => [t.id, i]))
  const rangChapitre = new Map(contenu.themes.flatMap((t) => t.chapitres.map((ch, i) => [ch.id, i] as const)))
  const dansFichier = new Map(contenu.cartes.map((c, i) => [c.id, i]))
  const parChapitre = [...contenu.cartes].sort(
    (a, b) =>
      (rangChapitre.get(a.chapitre) ?? 99) - (rangChapitre.get(b.chapitre) ?? 99) ||
      dansFichier.get(a.id)! - dansFichier.get(b.id)!,
  )
  const position = new Map<string, number>()
  const compteur = new Map<string, number>()
  for (const c of parChapitre) {
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
  if (f.etats.length && !f.etats.includes(etatAffiche(c, p, r))) return false
  return true
}

/** Notions apprises (ou introduites) aujourd'hui, pour le plafond quotidien. */
export function neuvesDuJour(contenu: Contenu, p: Progression, jour: string): number {
  return contenu.cartes.filter((c) => introduction(c, p) === jour).length
}

/** Cartes à revoir, les plus en retard d'abord. */
export function cartesARevoir(contenu: Contenu, p: Progression, r: ReglagesSeance, maintenant: number, jour: string): ElementSeance[] {
  const retard = (c: Carte, volets: Volet[]) =>
    Math.max(...volets.map((v) => {
      const e = p.volets[`${c.id}#${v}`]
      return e ? joursEntre(e.echeance, jour) : 0
    }))
  return contenu.cartes
    .filter((c) => carteRetenue(c, p, r))
    .map((c) => ({ c, volets: voletsDus(c, p, maintenant, jour) }))
    .filter((x) => x.volets.length > 0)
    .map((x) => ({ ...x, retard: retard(x.c, x.volets) }))
    .sort((a, b) => b.retard - a.retard)
    .map((x) => ({ carte: x.c.id, volets: x.volets, neuve: false }))
}

/**
 * Notions à apprendre, dans l'ordre d'introduction, au plus `max`. Une notion n'est proposée que si ses prérequis sont au
 * moins « en cours » (appris) ; une notion prise plus tôt dans la même liste compte comme apprise.
 */
export function aApprendre(contenu: Contenu, p: Progression, r: ReglagesSeance, max: number): Carte[] {
  const parId = new Map(contenu.cartes.map((c) => [c.id, c]))
  const choisies: Carte[] = []
  const prises = new Set<string>()
  const candidates = ordreIntroduction(contenu).filter((c) => etatCarte(c, p) === 'neuve' && carteRetenue(c, p, r))
  // Plusieurs passes : une carte débloquée par une carte choisie plus loin dans l'ordre peut entrer à la passe suivante.
  let ajout = true
  while (ajout && choisies.length < max) {
    ajout = false
    for (const c of candidates) {
      if (choisies.length >= max) break
      if (prises.has(c.id) || prerequisManquants(c, parId, p, r, prises).length) continue
      choisies.push(c)
      prises.add(c.id)
      ajout = true
    }
  }
  return choisies
}

export function plafondRestant(contenu: Contenu, p: Progression, r: ReglagesSeance, jour: string): number {
  return Math.max(0, r.plafondNeuves - neuvesDuJour(contenu, p, jour))
}

export interface Lecon {
  chapitre: Chapitre
  /** Notions de la leçon, dans l'ordre. */
  cartes: Carte[]
}

/**
 * Leçon du jour : les prochaines notions à apprendre (plafond quotidien), celles d'un même chapitre (le premier qui en
 * a) ; les autres viendront dans une leçon suivante. null s'il n'y a rien à apprendre aujourd'hui.
 */
export function leconDuJour(contenu: Contenu, p: Progression, r: ReglagesSeance, jour: string): Lecon | null {
  const notions = aApprendre(contenu, p, r, plafondRestant(contenu, p, r, jour))
  if (!notions.length) return null
  const chapitre = contenu.themes.flatMap((t) => t.chapitres).find((ch) => ch.id === notions[0].chapitre)
  if (!chapitre) return null
  return { chapitre, cartes: notions.filter((c) => c.chapitre === chapitre.id) }
}

export interface ResumeJour {
  aRevoir: number
  /** Notions qu'on peut encore apprendre aujourd'hui (plafond compris). */
  aApprendre: number
  lecon: Lecon | null
  /** Durée estimée de la séance de révision. */
  minutesRevision: number
  /** Durée estimée de la leçon du jour. */
  minutesLecon: number
}

export function resumeDuJour(
  contenu: Contenu,
  p: Progression,
  r: ReglagesSeance & Pick<Reglages, 'seanceMinutes'>,
  maintenant: number,
  jour: string,
): ResumeJour {
  const revisions = cartesARevoir(contenu, p, r, maintenant, jour)
  const lecon = leconDuJour(contenu, p, r, jour)
  const seance = Math.min(revisions.length, r.seanceCartes)
  return {
    aRevoir: revisions.length,
    aApprendre: aApprendre(contenu, p, r, plafondRestant(contenu, p, r, jour)).length,
    lecon,
    minutesRevision: Math.min(r.seanceMinutes, Math.ceil(seance * MINUTES_REVISION)),
    minutesLecon: Math.ceil((lecon?.cartes.length ?? 0) * MINUTES_NOTION),
  }
}

/** Séance de révision courte : au plus `seanceCartes` cartes déjà apprises, les plus en retard d'abord. */
export function construireSeance(contenu: Contenu, p: Progression, reglages: ReglagesSeance, maintenant: number, jour: string): ElementSeance[] {
  return cartesARevoir(contenu, p, reglages, maintenant, jour).slice(0, reglages.seanceCartes)
}
