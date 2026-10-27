// Types du contenu (content/) : voir CLAUDE.md, « Format du contenu ».

export const NIVEAUX = ['terminale', 'l1', 'l2', 'l3', 'agreg'] as const
export type Niveau = (typeof NIVEAUX)[number]

export const NOMS_NIVEAUX: Record<Niveau, string> = {
  terminale: 'Terminale',
  l1: 'L1',
  l2: 'L2',
  l3: 'L3',
  agreg: 'Agrégation',
}

/** Rang d'un niveau, du plancher (0, Terminale) au plafond (4, agrégation). */
export const rangNiveau = (n: Niveau): number => NIVEAUX.indexOf(n)

export type TypeCarte = 'theoreme' | 'definition' | 'methode'

export interface Trou {
  reponse: string
  leurres: string[]
}

export interface Affirmation {
  texte: string
  vrai: boolean
  contre_exemple?: string
  justification?: string
}

export interface Carte {
  id: string
  theme: string
  niveau: Niveau
  chapitre: string
  type: TypeCarte
  titre: string
  recto: string
  /** Texte de la leçon qui précède l'énoncé quand on découvre la notion (intuition, usage, méthode, piège). */
  cours?: string
  hypotheses?: string[]
  hypotheses_pieges?: string[]
  enonce: string
  enonce_trous?: string
  trous?: Trou[]
  demo?: string
  exemple?: string
  contre_exemple?: string
  affirmations?: Affirmation[]
  prerequis: string[]
  doute?: string
}

export interface ReponseFinale {
  question: string
  propositions: string[]
  bonne: number
}

/** Variante d'un exercice : mêmes question et difficulté, d'autres valeurs (énoncé, réponse, correction propres). */
export interface Variante {
  enonce: string
  indices?: string[]
  reponse_finale?: ReponseFinale
  correction: string
}

export interface Exercice {
  id: string
  carte: string
  difficulte: number
  titre: string
  enonce: string
  duree_min: number
  indices: string[]
  reponse_finale?: ReponseFinale
  correction: string
  /** Autres jeux de valeurs (la variante 0 est l'exercice lui-même). */
  variantes?: Variante[]
  doute?: string
}

/** Variante `k` d'un exercice (0 = l'exercice tel qu'écrit), complétée par les champs communs. */
export function variante(e: Exercice, k: number): Exercice {
  const v = k > 0 ? e.variantes?.[k - 1] : undefined
  if (!v) return e
  return { ...e, enonce: v.enonce, correction: v.correction, indices: v.indices ?? e.indices, reponse_finale: v.reponse_finale }
}

export const nombreVariantes = (e: Exercice): number => 1 + (e.variantes?.length ?? 0)

export interface Chapitre {
  id: string
  niveau: Niveau
  titre: string
  /** Notions à couvrir, une carte chacune : le plan de rédaction et la mesure de la couverture. */
  notions: string[]
  /** Nombre de cartes visées pour ce chapitre (= nombre de notions ; total du tableau de bord). */
  prevues: number
}

export interface Theme {
  id: string
  titre: string
  chapitres: Chapitre[]
}

/** Volets d'une carte notés en répétition espacée (les exercices n'y entrent pas). */
export const VOLETS = ['enonce', 'demo', 'contre'] as const
export type Volet = (typeof VOLETS)[number]

export const NOMS_VOLETS: Record<Volet, string> = { enonce: 'Énoncé', demo: 'Démonstration', contre: 'Contre-exemple' }

/** Volets présents sur une carte : l'énoncé toujours, la démonstration et le contre-exemple s'ils existent. */
export function voletsDe(c: Carte): Volet[] {
  const v: Volet[] = ['enonce']
  if (c.demo) v.push('demo')
  if (c.contre_exemple) v.push('contre')
  return v
}

/** Identifiant de l'état de répétition d'un volet : `idCarte#enonce`, `#demo`, `#contre`. */
export const cleVolet = (carte: string, volet: Volet): string => `${carte}#${volet}`

/** Introduction d'un chapitre, lue en tête de sa leçon (content/lecons/<chapitre>.json). */
export interface IntroLecon {
  chapitre: string
  introduction: string
  doute?: string
}

export interface Contenu {
  themes: Theme[]
  cartes: Carte[]
  exercices: Exercice[]
  /** Introductions de leçon, par identifiant de chapitre (facultatives). */
  lecons?: IntroLecon[]
}
