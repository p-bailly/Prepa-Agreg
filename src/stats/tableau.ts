// Tableau de bord : série de jours, cartes acquises, progression et niveau atteint, cartes revues par jour, réussite par
// thème, exercices réussis par difficulté, compte à rebours. Fonctions pures (jour passé en argument).

import { cleVolet, NIVEAUX, rangNiveau, voletsDe, type Contenu, type Niveau } from '../contenu/types'
import { etatExercice } from '../exercices/niveaux'
import type { Progression } from '../lib/progression'
import type { Reglages } from '../lib/reglages'
import { ajouterJours, joursEntre } from '../lib/temps'
import { etatAffiche, sousLeDepart } from '../repetition/etat'

type R = Pick<Reglages, 'niveauDepart' | 'niveauxInferieurs'>

const actif = (p: Progression, jour: string) => {
  const a = p.jours[jour]
  return !!a && a.cartes + a.volets + a.exercices > 0
}

/** Jours consécutifs d'activité finissant aujourd'hui (ou hier, si rien n'est encore fait aujourd'hui). */
export function serieJours(p: Progression, jour: string): number {
  let j = actif(p, jour) ? jour : ajouterJours(jour, -1)
  let n = 0
  while (actif(p, j)) {
    n++
    j = ajouterJours(j, -1)
  }
  return n
}

export function cartesAcquises(c: Contenu, p: Progression, r: R): number {
  return c.cartes.filter((k) => !(sousLeDepart(k, r) && r.niveauxInferieurs === 'masques') && etatAffiche(k, p, r) === 'acquise')
    .length
}

export interface ProgresNiveau {
  niveau: Niveau
  acquises: number
  /** Cartes rédigées à ce niveau. */
  existantes: number
  /** Cartes visées (themes.json), au moins les cartes existantes. */
  visees: number
  /** Part acquise, de 0 à 1, rapportée aux cartes visées. */
  part: number
  /** Niveau sous le niveau de départ. */
  statut: 'normal' | 'acquis' | 'masque'
}

export function progresNiveaux(c: Contenu, p: Progression, r: R): ProgresNiveau[] {
  return NIVEAUX.map((niveau) => {
    const cartes = c.cartes.filter((k) => k.niveau === niveau)
    const prevues = c.themes.flatMap((t) => t.chapitres).filter((ch) => ch.niveau === niveau).reduce((s, ch) => s + ch.prevues, 0)
    const visees = Math.max(prevues, cartes.length)
    const sous = rangNiveau(niveau) < rangNiveau(r.niveauDepart)
    if (sous) {
      const statut = r.niveauxInferieurs === 'acquis' ? 'acquis' : 'masque'
      return { niveau, acquises: cartes.length, existantes: cartes.length, visees, part: statut === 'acquis' ? 1 : 0, statut }
    }
    const acquises = cartes.filter((k) => etatAffiche(k, p, r) === 'acquise').length
    return { niveau, acquises, existantes: cartes.length, visees, part: visees ? acquises / visees : 0, statut: 'normal' }
  })
}

/** Le plus haut niveau dont au moins 80 % des cartes visées sont acquises ; null si aucun. */
export function niveauAtteint(niveaux: ProgresNiveau[]): Niveau | null {
  let atteint: Niveau | null = null
  for (const n of niveaux) if (n.statut !== 'masque' && n.visees > 0 && n.part >= 0.8) atteint = n.niveau
  return atteint
}

/** Cartes revues par jour, des `n` derniers jours (le dernier est aujourd'hui). */
export function cartesParJour(p: Progression, jour: string, n = 30): { jour: string; cartes: number }[] {
  return Array.from({ length: n }, (_, i) => {
    const j = ajouterJours(jour, i - n + 1)
    return { jour: j, cartes: p.jours[j]?.cartes ?? 0 }
  })
}

export interface ReussiteTheme {
  theme: string
  titre: string
  revisions: number
  reussites: number
  part: number
}

/** Part des révisions notées Bien ou Facile, par thème, les thèmes faibles en premier (thèmes jamais révisés exclus). */
export function reussiteParTheme(c: Contenu, p: Progression): ReussiteTheme[] {
  return c.themes
    .map((t) => {
      let revisions = 0
      let reussites = 0
      for (const k of c.cartes.filter((x) => x.theme === t.id)) {
        for (const v of voletsDe(k)) {
          const e = p.volets[cleVolet(k.id, v)]
          revisions += e?.revisions ?? 0
          reussites += e?.reussites ?? 0
        }
      }
      return { theme: t.id, titre: t.titre, revisions, reussites, part: revisions ? reussites / revisions : 0 }
    })
    .filter((x) => x.revisions > 0)
    .sort((a, b) => a.part - b.part || b.revisions - a.revisions)
}

/** Exercices réussis (au moins une fois) par difficulté, sur l'ensemble des exercices existants. */
export function exercicesParDifficulte(c: Contenu, p: Progression): { difficulte: number; reussis: number; total: number }[] {
  return [1, 2, 3, 4, 5].map((d) => {
    const ex = c.exercices.filter((e) => e.difficulte === d)
    return { difficulte: d, total: ex.length, reussis: ex.filter((e) => etatExercice(p.exercices[e.id]) === 'reussi').length }
  })
}

/** Jours restants avant l'écrit ; null sans date (la date est facultative et rien d'autre n'en dépend). */
export function joursAvantEcrit(dateEcrit: string | null, jour: string): number | null {
  return dateEcrit ? joursEntre(jour, dateEcrit) : null
}
