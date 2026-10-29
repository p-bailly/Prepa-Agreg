// Plan des chapitres : numérotation (dans l'ordre de l'accueil : niveau, puis thème, puis ordre de themes.json) et
// avancement d'une notion ou d'un chapitre, d'après les exercices réussis, le cours lu et le meilleur QCM.

import type { Exercice } from '../contenu/types'
import { NIVEAUX, type Carte, type Chapitre, type Contenu, type Theme } from '../contenu/types'
import { etatExercice } from '../exercices/niveaux'
import type { Progression } from '../lib/progression'

export interface ChapitreNumerote {
  numero: number
  chapitre: Chapitre
  theme: Theme
}

/** Tous les chapitres, numérotés de 1 à n : Terminale d'abord, puis L1…, et dans un niveau l'ordre des thèmes. */
export function chapitresNumerotes(c: Contenu): ChapitreNumerote[] {
  const sortie: ChapitreNumerote[] = []
  for (const n of NIVEAUX)
    for (const theme of c.themes)
      for (const chapitre of theme.chapitres) if (chapitre.niveau === n) sortie.push({ numero: sortie.length + 1, chapitre, theme })
  return sortie
}

/** Notions d'un chapitre (ses cartes rédigées), dans l'ordre des fichiers. */
export const notionsDe = (c: Contenu, chapitre: string): Carte[] => c.cartes.filter((k) => k.chapitre === chapitre)

/** À faire : rien d'ouvert ; en cours : cours lu ou exercice tenté ; fait : tous les exercices réussis (ou cours lu sans exercice). */
export type Avancement = 'vide' | 'a-faire' | 'en-cours' | 'fait'

export function avancementNotion(k: Carte, exercices: Exercice[], p: Progression): Exclude<Avancement, 'vide'> {
  const etats = exercices.map((e) => etatExercice(p.exercices[e.id]))
  if (etats.length ? etats.every((e) => e === 'reussi') : !!p.apprises[k.id]) return 'fait'
  return p.apprises[k.id] || etats.some((e) => e !== 'a-faire') ? 'en-cours' : 'a-faire'
}

/** Meilleur résultat au QCM du chapitre, ou null s'il n'a jamais été fait. */
export function meilleurQcm(p: Progression, chapitre: string): { bonnes: number; total: number } | null {
  const r = p.qcm[chapitre] ?? []
  return r.length ? r.reduce((a, b) => (b.bonnes / b.total > a.bonnes / a.total ? b : a)) : null
}

/** Seuil du QCM pour qu'un chapitre compte comme fait. */
export const SEUIL_QCM = 0.8

/** Vide : pas encore rédigé ; fait : toutes ses notions faites et au moins 80 % à un QCM. */
export function avancementChapitre(c: Contenu, chapitre: string, p: Progression): Avancement {
  const notions = notionsDe(c, chapitre)
  if (!notions.length) return 'vide'
  const etats = notions.map((k) => avancementNotion(k, c.exercices.filter((e) => e.carte === k.id), p))
  const qcm = meilleurQcm(p, chapitre)
  if (etats.every((e) => e === 'fait') && qcm && qcm.bonnes / qcm.total >= SEUIL_QCM) return 'fait'
  return qcm || etats.some((e) => e !== 'a-faire') ? 'en-cours' : 'a-faire'
}
