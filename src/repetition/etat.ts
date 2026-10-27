// État d'une carte (neuve, en cours, acquise), prérequis et déblocage.

import type { EtatCarte, Reglages } from '../lib/reglages'
import { cleVolet, rangNiveau, voletsDe, type Carte, type Volet } from '../contenu/types'
import { estDu, SEUIL_ACQUIS, type EtatVolet } from './planification'

type Volets = Record<string, EtatVolet>

/**
 * Neuve : aucun volet jamais noté. Acquise : chaque volet présent a un intervalle d'au moins trois semaines.
 * En cours : entre les deux.
 */
export function etatCarte(c: Carte, volets: Volets): EtatCarte {
  const etats = voletsDe(c).map((v) => volets[cleVolet(c.id, v)])
  if (etats.every((e) => !e)) return 'neuve'
  return etats.every((e) => e && e.intervalle >= SEUIL_ACQUIS) ? 'acquise' : 'en-cours'
}

/** Carte d'un niveau inférieur au niveau de départ (considérée comme acquise ou masquée selon les réglages). */
export const sousLeDepart = (c: Carte, r: Pick<Reglages, 'niveauDepart'>): boolean =>
  rangNiveau(c.niveau) < rangNiveau(r.niveauDepart)

/** État affiché : une carte sous le niveau de départ « considérée comme acquise » compte comme acquise. */
export function etatAffiche(c: Carte, volets: Volets, r: Pick<Reglages, 'niveauDepart' | 'niveauxInferieurs'>): EtatCarte {
  if (sousLeDepart(c, r) && r.niveauxInferieurs === 'acquis' && etatCarte(c, volets) === 'neuve') return 'acquise'
  return etatCarte(c, volets)
}

/**
 * Prérequis pas encore satisfaits. Un prérequis est satisfait dès qu'il est au moins « en cours », ou s'il est sous le
 * niveau de départ (acquis ou masqué : on ne bloque pas sur ce qu'on a choisi de ne pas réviser). `enCoursEnPlus` :
 * cartes qu'on s'apprête à introduire dans la même séance. Un prérequis inconnu (identifiant absent) bloque.
 */
export function prerequisManquants(
  c: Carte,
  parId: Map<string, Carte>,
  volets: Volets,
  r: Pick<Reglages, 'niveauDepart'>,
  enCoursEnPlus: ReadonlySet<string> = new Set(),
): string[] {
  return c.prerequis.filter((id) => {
    const p = parId.get(id)
    if (!p) return true
    if (sousLeDepart(p, r) || enCoursEnPlus.has(id)) return false
    return etatCarte(p, volets) === 'neuve'
  })
}

/**
 * Volets à revoir d'une carte déjà introduite, dans l'ordre énoncé, démonstration, contre-exemple. Un volet jamais noté
 * d'une carte introduite est à voir (la démonstration n'a pas forcément été ouverte le premier jour).
 */
export function voletsDus(c: Carte, volets: Volets, maintenant: number, jour: string): Volet[] {
  if (etatCarte(c, volets) === 'neuve') return []
  return voletsDe(c).filter((v) => {
    const e = volets[cleVolet(c.id, v)]
    return !e || estDu(e, maintenant, jour)
  })
}
