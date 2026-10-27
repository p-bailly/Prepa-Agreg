// Contenu factice pour les tests de la logique (indépendant de content/).
import type { Carte, Contenu, Exercice, Niveau } from '../contenu/types'
import { REGLAGES_DEFAUT, type Reglages } from '../lib/reglages'

export function carte(id: string, theme: string, niveau: Niveau, prerequis: string[] = [], extra: Partial<Carte> = {}): Carte {
  return {
    id,
    theme,
    niveau,
    chapitre: `${theme}-${niveau}`,
    type: 'theoreme',
    titre: id,
    recto: `Énoncer ${id}.`,
    enonce: `Énoncé de ${id}.`,
    prerequis,
    ...extra,
  }
}

export function exercice(id: string, carteId: string, difficulte: number): Exercice {
  return { id, carte: carteId, difficulte, titre: id, enonce: '…', duree_min: 5, indices: [], correction: '…' }
}

export function contenu(cartes: Carte[], exercices: Exercice[] = []): Contenu {
  const themes = [...new Set(cartes.map((c) => c.theme))].map((t) => ({
    id: t,
    titre: t,
    chapitres: [...new Set(cartes.filter((c) => c.theme === t).map((c) => c.niveau))].map((n) => ({
      id: `${t}-${n}`,
      niveau: n,
      titre: `${t} ${n}`,
      prevues: cartes.filter((c) => c.theme === t && c.niveau === n).length,
    })),
  }))
  return { themes, cartes, exercices }
}

export const reglages = (r: Partial<Reglages> = {}): Reglages => ({ ...REGLAGES_DEFAUT, ...r })
