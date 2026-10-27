// Chargement du contenu : tous les fichiers de content/ sont embarqués dans l'APK à la construction (import.meta.glob).
// Ajouter du contenu = ajouter un fichier content/cartes/<theme>/<niveau>.json ou content/exercices/<theme>/<niveau>.json.

import themes from '../../content/themes.json'
import type { Carte, Contenu, Exercice, IntroLecon, Theme } from './types'

const fichiersCartes = import.meta.glob<{ cartes: Carte[] }>('../../content/cartes/*/*.json', { eager: true, import: 'default' })
const fichiersExercices = import.meta.glob<{ exercices: Exercice[] }>('../../content/exercices/*/*.json', {
  eager: true,
  import: 'default',
})

const fichiersLecons = import.meta.glob<IntroLecon>('../../content/lecons/*.json', { eager: true, import: 'default' })

export const CONTENU: Contenu = {
  themes: (themes as { themes: Theme[] }).themes,
  cartes: Object.values(fichiersCartes).flatMap((f) => f.cartes),
  exercices: Object.values(fichiersExercices).flatMap((f) => f.exercices),
  lecons: Object.values(fichiersLecons),
}
