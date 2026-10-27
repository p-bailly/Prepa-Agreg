// Répétition espacée par volet, SM-2 simplifié. Chaque volet présent d'une carte (énoncé, démonstration, contre-exemple)
// a son propre état, sous la clé `idCarte#volet`. Quatre notes : À revoir, Difficile, Bien, Facile.
// Les échéances sont des jours civils (« AAAA-MM-JJ ») ; seule la reprise après « À revoir » est en minutes.

import { ajouterJours, joursEntre } from '../lib/temps'

export type Note = 0 | 1 | 2 | 3
export const NOTES: Note[] = [0, 1, 2, 3]
export const NOMS_NOTES: Record<Note, string> = { 0: 'À revoir', 1: 'Difficile', 2: 'Bien', 3: 'Facile' }

export interface EtatVolet {
  /** Intervalle en jours depuis la dernière révision (0 : en apprentissage, revu dans la séance). */
  intervalle: number
  /** Facteur de facilité SM-2, entre FACILITE_MIN et FACILITE_MAX. */
  facilite: number
  /** Révisions réussies d'affilée (Difficile ou mieux). */
  repetitions: number
  /** Nombre de fois où le volet a été oublié (« À revoir » après une première réussite). */
  echecs: number
  /** Jour à partir duquel le volet est à revoir. */
  echeance: string
  /** Instant de reprise après « À revoir » (dans la séance ou plus tard le même jour). */
  echeanceMs?: number
  /** Instant de la dernière révision. */
  derniere: number
  /** Jour de la première révision (sert au plafond quotidien de cartes neuves). */
  premiere: string
  revisions: number
  /** Révisions notées Bien ou Facile. */
  reussites: number
}

export const FACILITE_DEPART = 2.5
export const FACILITE_MIN = 1.3
export const FACILITE_MAX = 3
export const INTERVALLE_MAX = 365
/** Délai de reprise après « À revoir ». */
export const REPRISE_MINUTES = 10
/** Un volet est acquis quand son intervalle atteint trois semaines. */
export const SEUIL_ACQUIS = 21

/** Intervalles (jours) pour Difficile, Bien, Facile la toute première fois qu'on voit un volet. */
const PREMIERE_FOIS: Record<1 | 2 | 3, number> = { 1: 1, 2: 4, 3: 9 }
/** Même chose pour un volet oublié puis repris (on repart plus prudemment). */
const REPRISE: Record<1 | 2 | 3, number> = { 1: 1, 2: 2, 3: 4 }

const borner = (x: number, min: number, max: number) => Math.min(max, Math.max(min, x))

/**
 * Nouvel état d'un volet après une note. `maintenant` est l'instant de la révision et `jour` le jour civil local
 * correspondant (calculé par l'appelant avec l'horloge injectée).
 */
export function planifier(etat: EtatVolet | undefined, note: Note, maintenant: number, jour: string): EtatVolet {
  const avant: EtatVolet = etat ?? {
    intervalle: 0,
    facilite: FACILITE_DEPART,
    repetitions: 0,
    echecs: 0,
    echeance: jour,
    derniere: maintenant,
    premiere: jour,
    revisions: 0,
    reussites: 0,
  }
  const commun = {
    derniere: maintenant,
    premiere: avant.premiere,
    revisions: avant.revisions + 1,
    reussites: avant.reussites + (note >= 2 ? 1 : 0),
  }
  if (note === 0) {
    return {
      ...commun,
      intervalle: 0,
      facilite: borner(avant.facilite - 0.2, FACILITE_MIN, FACILITE_MAX),
      repetitions: 0,
      // Oublier un volet déjà réussi compte comme un échec ; rater une carte qu'on découvre, non.
      echecs: avant.echecs + (avant.intervalle > 0 ? 1 : 0),
      echeance: jour,
      echeanceMs: maintenant + REPRISE_MINUTES * 60_000,
    }
  }
  const facilite = borner(avant.facilite + (note === 1 ? -0.15 : note === 3 ? 0.15 : 0), FACILITE_MIN, FACILITE_MAX)
  const intervalle = Math.min(INTERVALLE_MAX, prochainIntervalle(avant, note, jour))
  return {
    ...commun,
    intervalle,
    facilite,
    repetitions: avant.repetitions + 1,
    echecs: avant.echecs,
    echeance: ajouterJours(jour, intervalle),
  }
}

function prochainIntervalle(e: EtatVolet, note: 1 | 2 | 3, jour: string): number {
  if (e.intervalle === 0) return (e.revisions === 0 || e.echecs === 0 ? PREMIERE_FOIS : REPRISE)[note]
  // Une carte revue en retard a tenu plus longtemps que prévu : on en tient compte (à moitié pour « Bien »).
  const retard = Math.max(0, joursEntre(e.echeance, jour))
  const I = e.intervalle
  const difficile = Math.max(I + 1, Math.round(I * 1.2))
  if (note === 1) return difficile
  const bien = Math.max(difficile + 1, Math.round((I + retard / 2) * e.facilite))
  if (note === 2) return bien
  return Math.max(bien + 1, Math.round((I + retard) * e.facilite * 1.3))
}

/** Le volet est-il à revoir à cet instant ? (Un volet jamais vu n'a pas d'état : c'est à l'appelant d'en décider.) */
export function estDu(e: EtatVolet, maintenant: number, jour: string): boolean {
  if (e.echeanceMs !== undefined) return maintenant >= e.echeanceMs || joursEntre(e.echeance, jour) > 0
  return joursEntre(e.echeance, jour) >= 0
}

/** Délai affiché sous un bouton de notation : « < 10 min », « 1 j », « 4 j », « 2 mois », « 1 an ». */
export function delaiAffiche(etat: EtatVolet | undefined, note: Note, maintenant: number, jour: string): string {
  if (note === 0) return `< ${REPRISE_MINUTES} min`
  return formaterJours(planifier(etat, note, maintenant, jour).intervalle)
}

export function formaterJours(j: number): string {
  if (j < 30) return `${j} j`
  if (j < 365) return `${Math.round(j / 30)} mois`
  return '1 an'
}
