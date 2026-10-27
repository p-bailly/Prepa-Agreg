// Export et import de la progression en JSON : le stockage est local, c'est la seule sauvegarde. L'export contient aussi
// les réglages et les signalements d'erreurs (dans la progression).
import { normaliserProgression, type Progression } from './progression'
import { normaliserReglages, type Reglages } from './reglages'

export const FORMAT_SAUVEGARDE = 'prepa-agreg/progression-v1'

export interface Sauvegarde {
  format: typeof FORMAT_SAUVEGARDE
  exporte: string
  version: string
  reglages: Reglages
  progression: Progression
}

export function fabriquerSauvegarde(p: Progression, r: Reglages, maintenant: number, version: string): Sauvegarde {
  return { format: FORMAT_SAUVEGARDE, exporte: new Date(maintenant).toISOString(), version, reglages: r, progression: p }
}

/** Lit un fichier de sauvegarde ; lève une erreur au message lisible s'il n'en est pas un. */
export function lireSauvegarde(texte: string): { progression: Progression; reglages: Reglages } {
  let brut: unknown
  try {
    brut = JSON.parse(texte)
  } catch {
    throw new Error('Ce fichier n’est pas du JSON.')
  }
  const s = brut as Partial<Sauvegarde> | null
  if (!s || typeof s !== 'object' || s.format !== FORMAT_SAUVEGARDE)
    throw new Error('Ce fichier n’est pas une sauvegarde de Prépa Agreg (format progression-v1).')
  if (!s.progression || typeof s.progression !== 'object') throw new Error('La sauvegarde ne contient pas de progression.')
  return { progression: normaliserProgression(s.progression), reglages: normaliserReglages(s.reglages) }
}
