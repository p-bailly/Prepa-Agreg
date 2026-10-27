// Les quatre modes de révision de l'énoncé, choisis selon les données de la carte et alternés d'une carte à l'autre.

import type { Carte } from '../contenu/types'
import { MODES, type Mode } from '../lib/reglages'

/** Modes qu'une carte permet : un mode n'apparaît que si la carte a ce qu'il faut. */
export function modesPossibles(c: Carte): Mode[] {
  const m: Mode[] = ['tete']
  if (c.hypotheses?.length && c.hypotheses_pieges?.length) m.push('hypotheses')
  if (c.enonce_trous && c.trous?.length) m.push('trous')
  if (c.affirmations?.some((a) => !a.vrai && a.contre_exemple)) m.push('vraifaux')
  return m
}

/**
 * Mode de la prochaine carte : le premier mode possible et actif qui suit `dernier` dans l'ordre des modes (on alterne).
 * Une carte neuve se découvre toujours en « énoncer de tête » ; si aucun mode actif ne convient, on revient à celui-là.
 */
export function choisirMode(c: Carte, actifs: readonly Mode[], dernier: Mode | null, neuve: boolean): Mode {
  if (neuve) return 'tete'
  const possibles = modesPossibles(c).filter((m) => actifs.includes(m))
  if (!possibles.length) return 'tete'
  const depart = dernier ? MODES.indexOf(dernier) : -1
  for (let k = 1; k <= MODES.length; k++) {
    const m = MODES[(depart + k) % MODES.length]
    if (possibles.includes(m)) return m
  }
  return possibles[0]
}
