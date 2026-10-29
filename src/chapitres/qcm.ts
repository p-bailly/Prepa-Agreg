// QCM d'un chapitre, construit à partir du contenu déjà vérifié (voir CLAUDE.md, « Contrôle de l'exactitude ») :
// affirmations des cartes (vrai ou faux, avec justification ou contre-exemple), énoncés à trous (le bon morceau parmi
// quatre, les leurres rendant l'énoncé faux) et questions des exercices de difficulté 1 ou 2 (avec leur énoncé).

import type { Contenu } from '../contenu/types'
import { melanger } from '../lib/hasard'

export type GenreQuestion = 'vrai-faux' | 'completer' | 'calcul'

export interface QuestionQcm {
  id: string
  genre: GenreQuestion
  /** Notion (carte) dont vient la question. */
  carte: string
  /** Énoncé de l'exercice, pour une question de calcul. */
  contexte?: string
  enonce: string
  propositions: string[]
  bonne: number
  /** Justification, contre-exemple ou correction, montrés après la réponse. */
  explication?: string
}

export const CONSIGNES: Record<GenreQuestion, string> = {
  'vrai-faux': 'Vrai ou faux ?',
  completer: 'Compléter l’énoncé',
  calcul: 'Question d’exercice',
}

/** Toutes les questions possibles d'un chapitre. */
export function questionsChapitre(c: Contenu, chapitre: string): QuestionQcm[] {
  const sortie: QuestionQcm[] = []
  const cartes = c.cartes.filter((k) => k.chapitre === chapitre)
  for (const k of cartes) {
    k.affirmations?.forEach((a, i) =>
      sortie.push({
        id: `${k.id}#affirmation${i}`,
        genre: 'vrai-faux',
        carte: k.id,
        enonce: a.texte,
        propositions: ['Vrai', 'Faux'],
        bonne: a.vrai ? 0 : 1,
        explication: a.vrai ? a.justification : a.contre_exemple && `**Contre-exemple.** ${a.contre_exemple}`,
      }),
    )
    if (k.enonce_trous && k.trous)
      k.trous.forEach((t, i) => {
        const texte = k.enonce_trous!.replace(/\{\{(\d+)\}\}/g, (_, n: string) =>
          Number(n) === i + 1 ? '**[ … ]**' : (k.trous![Number(n) - 1]?.reponse ?? '…'),
        )
        sortie.push({
          id: `${k.id}#trou${i}`,
          genre: 'completer',
          carte: k.id,
          enonce: texte,
          propositions: [t.reponse, ...t.leurres],
          bonne: 0,
          explication: `**Énoncé complet.** ${k.enonce}`,
        })
      })
  }
  const ids = new Set(cartes.map((k) => k.id))
  for (const e of c.exercices)
    if (ids.has(e.carte) && e.difficulte <= 2 && e.reponse_finale)
      sortie.push({
        id: `${e.id}#reponse`,
        genre: 'calcul',
        carte: e.carte,
        contexte: e.enonce,
        enonce: e.reponse_finale.question,
        propositions: e.reponse_finale.propositions,
        bonne: e.reponse_finale.bonne,
        explication: e.correction,
      })
  return sortie
}

/**
 * Tire `n` questions en alternant les genres (pour ne pas avoir que des vrai-faux), dans un ordre mélangé, propositions
 * mélangées sauf « Vrai / Faux » qui garde son ordre.
 */
export function tirerQcm(questions: QuestionQcm[], n = 10, alea: () => number = Math.random): QuestionQcm[] {
  const parGenre = (['vrai-faux', 'completer', 'calcul'] as GenreQuestion[]).map((g) => melanger(questions.filter((q) => q.genre === g), alea))
  const choisies: QuestionQcm[] = []
  for (let tour = 0; choisies.length < n && parGenre.some((l) => l.length > tour); tour++)
    for (const l of parGenre) if (l[tour] && choisies.length < n) choisies.push(l[tour])
  return melanger(choisies, alea).map((q) => {
    if (q.genre === 'vrai-faux') return q
    const ordre = melanger(q.propositions.map((_, i) => i), alea)
    return { ...q, propositions: ordre.map((i) => q.propositions[i]), bonne: ordre.indexOf(q.bonne) }
  })
}
