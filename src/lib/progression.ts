// Progression, sur l'appareil sous la clé versionnée `progression-v1` (changer le format = changer le suffixe).
// Indexée par identifiant de carte (volets `idCarte#volet`) et d'exercice : renommer un identifiant perd son historique,
// donc ne jamais le faire sans migration. Fonctions pures : l'instant et le jour sont passés en argument.

import { cleVolet, type Volet } from '../contenu/types'
import { planifier, type EtatVolet, type Note } from '../repetition/planification'
import type { Mode } from './reglages'
import { ecrire, lire } from './stockage'

export type ResultatExercice = 'reussi' | 'partiel' | 'rate'

export interface Tentative {
  /** Instant de l'auto-évaluation. */
  date: number
  jour: string
  resultat: ResultatExercice
  /** Indices débloqués avant de conclure. */
  indices: number
  /** Réponse choisie parmi les propositions, si l'exercice en a. */
  reponse?: 'juste' | 'fausse'
  /** Variante faite (0 = l'exercice tel qu'écrit ; absent dans les anciennes tentatives). */
  variante?: number
}

export interface Signalement {
  date: number
  /** Identifiant de la carte ou de l'exercice. */
  cible: string
  genre: 'carte' | 'exercice'
  /** Volet ou écran concerné (« énoncé », « correction »…), si connu. */
  ou?: string
  commentaire: string
}

export interface ActiviteJour {
  /** Cartes revues (une par carte présentée et notée). */
  cartes: number
  /** Notes données à des volets (toutes). */
  volets: number
  exercices: number
}

export interface Progression {
  volets: Record<string, EtatVolet>
  exercices: Record<string, Tentative[]>
  jours: Record<string, ActiviteJour>
  signalements: Signalement[]
  /** Dernier mode de révision utilisé, pour les alterner d'une carte à l'autre. */
  dernierMode: Mode | null
}

export const PROGRESSION_VIDE: Progression = { volets: {}, exercices: {}, jours: {}, signalements: [], dernierMode: null }
export const CLE_PROGRESSION = 'progression-v1'

export function normaliserProgression(brut: unknown): Progression {
  const p = { ...PROGRESSION_VIDE, ...(brut && typeof brut === 'object' ? (brut as Partial<Progression>) : {}) }
  const objet = <T,>(x: unknown): Record<string, T> =>
    x && typeof x === 'object' && !Array.isArray(x) ? (x as Record<string, T>) : {}
  return {
    volets: objet<EtatVolet>(p.volets),
    exercices: objet<Tentative[]>(p.exercices),
    jours: objet<ActiviteJour>(p.jours),
    signalements: Array.isArray(p.signalements) ? p.signalements : [],
    dernierMode: typeof p.dernierMode === 'string' ? p.dernierMode : null,
  }
}

export const lireProgression = (): Progression => normaliserProgression(lire(CLE_PROGRESSION, {}))
export const ecrireProgression = (p: Progression): void => ecrire(CLE_PROGRESSION, p)

function ajouterActivite(p: Progression, jour: string, delta: Partial<ActiviteJour>): Record<string, ActiviteJour> {
  const a = p.jours[jour] ?? { cartes: 0, volets: 0, exercices: 0 }
  return {
    ...p.jours,
    [jour]: {
      cartes: a.cartes + (delta.cartes ?? 0),
      volets: a.volets + (delta.volets ?? 0),
      exercices: a.exercices + (delta.exercices ?? 0),
    },
  }
}

/**
 * Note un volet. `compterCarte` : vrai pour la note qui conclut la présentation d'une carte (le volet principal),
 * pour que « cartes revues par jour » compte des cartes et non des volets.
 */
export function noterVolet(
  p: Progression,
  carte: string,
  volet: Volet,
  note: Note,
  maintenant: number,
  jour: string,
  compterCarte: boolean,
): Progression {
  const cle = cleVolet(carte, volet)
  return {
    ...p,
    volets: { ...p.volets, [cle]: planifier(p.volets[cle], note, maintenant, jour) },
    jours: ajouterActivite(p, jour, { volets: 1, cartes: compterCarte ? 1 : 0 }),
  }
}

export function enregistrerTentative(p: Progression, exercice: string, t: Tentative): Progression {
  return {
    ...p,
    exercices: { ...p.exercices, [exercice]: [...(p.exercices[exercice] ?? []), t] },
    jours: ajouterActivite(p, t.jour, { exercices: 1 }),
  }
}

export function ajouterSignalement(p: Progression, s: Signalement): Progression {
  return { ...p, signalements: [...p.signalements, s] }
}

export function retenirMode(p: Progression, mode: Mode): Progression {
  return p.dernierMode === mode ? p : { ...p, dernierMode: mode }
}
