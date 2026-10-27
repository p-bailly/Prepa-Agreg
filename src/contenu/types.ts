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
  doute?: string
}

export interface Chapitre {
  id: string
  niveau: Niveau
  titre: string
  /** Nombre de cartes visées pour ce chapitre (total du tableau de bord). */
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

export interface Contenu {
  themes: Theme[]
  cartes: Carte[]
  exercices: Exercice[]
}
