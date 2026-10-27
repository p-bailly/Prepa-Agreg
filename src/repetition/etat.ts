// État d'une carte (neuve, en cours, acquise), prérequis et déblocage.

import type { EtatCarte, Reglages } from '../lib/reglages'
import { cleVolet, rangNiveau, voletsDe, type Carte, type Volet } from '../contenu/types'
import { estDu, SEUIL_ACQUIS, type EtatVolet } from './planification'

/** Ce que l'état d'une carte regarde dans la progression. */
export interface Suivi {
  volets: Record<string, EtatVolet>
  apprises: Record<string, string>
}

/**
 * Neuve : ni apprise en leçon, ni jamais notée. Acquise : chaque volet présent a un intervalle d'au moins trois semaines.
 * En cours : entre les deux (apprise, ou en révision).
 */
export function etatCarte(c: Carte, suivi: Suivi): EtatCarte {
  const etats = voletsDe(c).map((v) => suivi.volets[cleVolet(c.id, v)])
  if (etats.every((e) => !e)) return suivi.apprises[c.id] ? 'en-cours' : 'neuve'
  return etats.every((e) => e && e.intervalle >= SEUIL_ACQUIS) ? 'acquise' : 'en-cours'
}

/** Carte d'un niveau inférieur au niveau de départ (considérée comme acquise ou masquée selon les réglages). */
export const sousLeDepart = (c: Carte, r: Pick<Reglages, 'niveauDepart'>): boolean =>
  rangNiveau(c.niveau) < rangNiveau(r.niveauDepart)

/** État affiché : une carte sous le niveau de départ « considérée comme acquise » compte comme acquise. */
export function etatAffiche(c: Carte, suivi: Suivi, r: Pick<Reglages, 'niveauDepart' | 'niveauxInferieurs'>): EtatCarte {
  if (sousLeDepart(c, r) && r.niveauxInferieurs === 'acquis' && etatCarte(c, suivi) === 'neuve') return 'acquise'
  return etatCarte(c, suivi)
}

/**
 * Prérequis pas encore satisfaits. Un prérequis est satisfait dès qu'il est au moins « en cours » (appris), ou s'il est
 * sous le niveau de départ (acquis ou masqué : on ne bloque pas sur ce qu'on a choisi de ne pas réviser). `enCoursEnPlus` :
 * cartes qu'on s'apprête à apprendre dans la même leçon. Un prérequis inconnu (identifiant absent) bloque.
 */
export function prerequisManquants(
  c: Carte,
  parId: Map<string, Carte>,
  suivi: Suivi,
  r: Pick<Reglages, 'niveauDepart'>,
  enCoursEnPlus: ReadonlySet<string> = new Set(),
): string[] {
  return c.prerequis.filter((id) => {
    const p = parId.get(id)
    if (!p) return true
    if (sousLeDepart(p, r) || enCoursEnPlus.has(id)) return false
    return etatCarte(p, suivi) === 'neuve'
  })
}

/** Jour où la carte a été apprise en leçon ou, à défaut, notée pour la première fois (null si neuve). */
export function introduction(c: Carte, suivi: Suivi): string | null {
  const jours = [suivi.apprises[c.id], ...voletsDe(c).map((v) => suivi.volets[cleVolet(c.id, v)]?.premiere)]
    .filter((x): x is string => !!x)
    .sort()
  return jours[0] ?? null
}

/**
 * Volets à revoir d'une carte apprise, dans l'ordre énoncé, démonstration, contre-exemple. Un volet jamais noté (carte
 * apprise en leçon, ou démonstration pas ouverte) est à voir à partir du lendemain de l'apprentissage : on ne révise
 * jamais le jour même ce qu'on vient d'apprendre.
 */
export function voletsDus(c: Carte, suivi: Suivi, maintenant: number, jour: string): Volet[] {
  const debut = introduction(c, suivi)
  if (debut === null) return []
  return voletsDe(c).filter((v) => {
    const e = suivi.volets[cleVolet(c.id, v)]
    return e ? estDu(e, maintenant, jour) : jour > debut
  })
}
