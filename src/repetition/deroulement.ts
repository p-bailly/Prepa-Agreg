// Déroulé d'une séance courte : avancer après la note du volet principal, remettre en fin de séance une carte notée
// « À revoir » (au plus deux fois), arrêter quand il n'y a plus de carte ou que la durée est écoulée.

import type { Volet } from '../contenu/types'
import type { Mode } from '../lib/reglages'
import type { Note } from './planification'
import type { ElementSeance } from './seance'

export interface EtatSeance {
  elements: ElementSeance[]
  index: number
  debut: number
  /** Notes des volets principaux, pour le bilan. */
  notes: Note[]
  /** Cartes distinctes vues. */
  vues: string[]
  /** Mode retenu pour la carte en cours (fixé à son arrivée). */
  mode: Mode | null
  fin: boolean
}

export const REPRISES_MAX = 2

export const nouvelleSeance = (elements: ElementSeance[], debut: number): EtatSeance => ({
  elements,
  index: 0,
  debut,
  notes: [],
  vues: [],
  mode: null,
  fin: elements.length === 0,
})

/** État après la note du volet principal de la carte en cours, à l'instant `t`. */
export function avancer(s: EtatSeance, volet: Volet, note: Note, t: number, dureeMinutes: number): EtatSeance {
  const element = s.elements[s.index]
  if (!element || s.fin) return s
  const carte = element.carte
  const vues = s.vues.includes(carte) ? s.vues : [...s.vues, carte]
  const passages = s.elements.slice(0, s.index + 1).filter((e) => e.carte === carte).length
  const elements =
    note === 0 && passages <= REPRISES_MAX ? [...s.elements, { carte, volets: [volet], neuve: false }] : s.elements
  const index = s.index + 1
  // La carte en cours est toujours terminée ; la séance s'arrête ensuite si le temps est écoulé.
  const tempsEcoule = t - s.debut >= dureeMinutes * 60_000
  return { ...s, elements, index, vues, notes: [...s.notes, note], mode: null, fin: index >= elements.length || tempsEcoule }
}

/** « n / N » : rang de la carte en cours parmi les cartes distinctes de la séance (les reprises ne comptent pas). */
export function position(s: EtatSeance): { n: number; total: number } {
  const total = new Set(s.elements.map((e) => e.carte)).size
  const carte = s.elements[s.index]?.carte
  const n = Math.min(total, s.vues.length + (carte && s.vues.includes(carte) ? 0 : 1))
  return { n, total }
}
