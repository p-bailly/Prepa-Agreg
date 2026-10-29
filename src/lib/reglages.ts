// Réglages, tous sur l'appareil, sous la clé versionnée `reglages-v1` (changer le format = changer le suffixe).

import { NIVEAUX, type Niveau } from '../contenu/types'
import { ecrire, lire } from './stockage'

export const MODES = ['tete', 'hypotheses', 'trous', 'vraifaux'] as const
export type Mode = (typeof MODES)[number]
export const NOMS_MODES: Record<Mode, string> = {
  tete: 'Énoncer de tête, puis révéler',
  hypotheses: 'Hypothèses à cocher',
  trous: 'Énoncé à trous',
  vraifaux: 'Vrai ou faux avec contre-exemple',
}

export const ETATS_CARTE = ['neuve', 'en-cours', 'acquise'] as const
export type EtatCarte = (typeof ETATS_CARTE)[number]
export const NOMS_ETATS: Record<EtatCarte, string> = { neuve: 'Neuves', 'en-cours': 'En cours', acquise: 'Acquises' }

export interface Filtres {
  /** Listes vides = pas de filtre. */
  themes: string[]
  niveaux: Niveau[]
  etats: EtatCarte[]
}

export interface Reglages {
  niveauDepart: Niveau
  /** Les niveaux sous le niveau de départ sont considérés comme acquis, ou masqués. */
  niveauxInferieurs: 'acquis' | 'masques'
  niveauxActifs: Niveau[]
  modesActifs: Mode[]
  plafondNeuves: number
  seanceCartes: number
  seanceMinutes: number
  /** Cartes à revoir par jour visées (affiché sur l'accueil). */
  objectifQuotidien: number
  /** Date de l'écrit « AAAA-MM-JJ » ; plus affichée ni modifiable depuis le 2026-10-29, gardée dans le format. */
  dateEcrit: string | null
  theme: 'auto' | 'clair' | 'sombre'
  /** Rappel quotidien retiré le 2026-10-29 : `actif` est remis à faux au démarrage (voir `App.tsx`). */
  rappel: { actif: boolean; heure: string }
  filtres: Filtres
}

export const REGLAGES_DEFAUT: Reglages = {
  niveauDepart: 'terminale',
  niveauxInferieurs: 'acquis',
  niveauxActifs: [...NIVEAUX],
  modesActifs: [...MODES],
  plafondNeuves: 10,
  seanceCartes: 15,
  seanceMinutes: 15,
  objectifQuotidien: 15,
  dateEcrit: null,
  theme: 'auto',
  rappel: { actif: false, heure: '19:00' },
  filtres: { themes: [], niveaux: [], etats: [] },
}

export const CLE_REGLAGES = 'reglages-v1'

/** Remet d'aplomb des réglages lus (fichier importé, ancienne version) : tout champ invalide reprend sa valeur par défaut. */
export function normaliserReglages(brut: unknown): Reglages {
  const r = { ...REGLAGES_DEFAUT, ...(brut && typeof brut === 'object' ? (brut as Partial<Reglages>) : {}) }
  const entier = (x: unknown, min: number, max: number, d: number) =>
    typeof x === 'number' && Number.isFinite(x) ? Math.min(max, Math.max(min, Math.round(x))) : d
  const liste = <T,>(x: unknown, permis: readonly T[], d: T[]): T[] =>
    Array.isArray(x) ? x.filter((v): v is T => permis.includes(v as T)) : d
  const niveauxActifs = liste(r.niveauxActifs, NIVEAUX, [...NIVEAUX])
  const modesActifs = liste(r.modesActifs, MODES, [...MODES])
  const f = r.filtres && typeof r.filtres === 'object' ? r.filtres : REGLAGES_DEFAUT.filtres
  return {
    niveauDepart: NIVEAUX.includes(r.niveauDepart) ? r.niveauDepart : 'terminale',
    niveauxInferieurs: r.niveauxInferieurs === 'masques' ? 'masques' : 'acquis',
    niveauxActifs: niveauxActifs.length ? niveauxActifs : [...NIVEAUX],
    modesActifs: modesActifs.length ? modesActifs : ['tete'],
    plafondNeuves: entier(r.plafondNeuves, 0, 100, 10),
    seanceCartes: entier(r.seanceCartes, 1, 100, 15),
    seanceMinutes: entier(r.seanceMinutes, 1, 120, 15),
    objectifQuotidien: entier(r.objectifQuotidien, 1, 500, 15),
    dateEcrit: typeof r.dateEcrit === 'string' && /^\d{4}-\d{2}-\d{2}$/.test(r.dateEcrit) ? r.dateEcrit : null,
    theme: r.theme === 'clair' || r.theme === 'sombre' ? r.theme : 'auto',
    rappel: {
      actif: !!r.rappel?.actif,
      heure: typeof r.rappel?.heure === 'string' && /^\d{2}:\d{2}$/.test(r.rappel.heure) ? r.rappel.heure : '19:00',
    },
    filtres: {
      themes: Array.isArray(f.themes) ? f.themes.filter((t): t is string => typeof t === 'string') : [],
      niveaux: liste(f.niveaux, NIVEAUX, []),
      etats: liste(f.etats, ETATS_CARTE, []),
    },
  }
}

export const lireReglages = (): Reglages => normaliserReglages(lire(CLE_REGLAGES, {}))
export const ecrireReglages = (r: Reglages): void => ecrire(CLE_REGLAGES, r)

/** Même calcul que le script d'index.html (thème posé avant le premier rendu). */
export function estSombre(r: Reglages): boolean {
  if (r.theme !== 'auto') return r.theme === 'sombre'
  return typeof matchMedia === 'function' && matchMedia('(prefers-color-scheme: dark)').matches
}
