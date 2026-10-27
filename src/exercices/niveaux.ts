// Exercices d'une notion : état des cinq niveaux de difficulté et exercice suggéré. Les exercices n'entrent pas dans la
// répétition espacée : on garde seulement l'historique des tentatives (date, résultat, indices utilisés).

import type { Exercice } from '../contenu/types'
import type { Tentative } from '../lib/progression'

/** ■ réussi, ◧ fait en partie, □ à faire, grisé : aucun exercice à ce niveau. */
export type EtatNiveau = 'reussi' | 'partiel' | 'a-faire' | 'aucun'
export const SYMBOLES: Record<EtatNiveau, string> = { reussi: '■', partiel: '◧', 'a-faire': '□', aucun: '□' }
export const NOMS_ETATS_NIVEAU: Record<EtatNiveau, string> = {
  reussi: 'réussi',
  partiel: 'fait en partie',
  'a-faire': 'à faire',
  aucun: 'aucun exercice',
}

export type EtatExercice = 'reussi' | 'partiel' | 'a-faire'

/** Un exercice est réussi s'il l'a été au moins une fois ; fait en partie s'il a été tenté sans réussite. */
export function etatExercice(historique: Tentative[] | undefined): EtatExercice {
  if (!historique?.length) return 'a-faire'
  return historique.some((t) => t.resultat === 'reussi') ? 'reussi' : 'partiel'
}

export interface Niveau {
  difficulte: number
  etat: EtatNiveau
  exercices: Exercice[]
}

/**
 * Les cinq niveaux d'une notion. Un niveau est réussi quand tous ses exercices ont été réussis au moins une fois, fait
 * en partie dès qu'un exercice y a été tenté, à faire sinon ; « aucun » s'il n'a pas d'exercice.
 */
export function niveauxDe(exercices: Exercice[], historique: Record<string, Tentative[]>): Niveau[] {
  return [1, 2, 3, 4, 5].map((d) => {
    const ex = exercices.filter((e) => e.difficulte === d)
    if (!ex.length) return { difficulte: d, etat: 'aucun', exercices: [] }
    const etats = ex.map((e) => etatExercice(historique[e.id]))
    const etat: EtatNiveau = etats.every((x) => x === 'reussi')
      ? 'reussi'
      : etats.some((x) => x !== 'a-faire')
        ? 'partiel'
        : 'a-faire'
    return { difficulte: d, etat, exercices: ex }
  })
}

/**
 * Variante à proposer : la première jamais faite ; si toutes l'ont été, celle faite il y a le plus longtemps.
 * Une tentative sans numéro de variante compte pour la variante 0.
 */
export function varianteSuivante(nombre: number, historique: Tentative[] | undefined): number {
  const derniere = new Array<number>(nombre).fill(-1)
  ;(historique ?? []).forEach((t, i) => {
    const k = t.variante ?? 0
    if (k >= 0 && k < nombre) derniere[k] = i
  })
  let choix = 0
  for (let k = 0; k < nombre; k++) if (derniere[k] < derniere[choix]) choix = k
  return choix
}

/** Premier niveau non réussi, et dans ce niveau le premier exercice non réussi ; null si tout est réussi. */
export function exerciceSuggere(exercices: Exercice[], historique: Record<string, Tentative[]>): Exercice | null {
  for (const n of niveauxDe(exercices, historique)) {
    if (n.etat === 'reussi' || n.etat === 'aucun') continue
    return n.exercices.find((e) => etatExercice(historique[e.id]) !== 'reussi') ?? null
  }
  return null
}
