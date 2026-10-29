// QCM d'un chapitre, construit à partir du contenu vérifié (voir CLAUDE.md, « Contrôle de l'exactitude ») : affirmations
// des cartes (vrai ou faux, avec justification ou contre-exemple), énoncés à trous (le bon morceau parmi quatre, les
// leurres rendant l'énoncé faux), questions des exercices de difficulté 1 ou 2 (avec leur énoncé) et questions écrites
// pour le chapitre (content/qcm/<chapitre>.json, vrai-faux ou choix parmi quatre). Un QCM porte sur un chapitre, ou sur
// tout un niveau (questions à choix d'un côté, vrai-faux de l'autre, réparties entre les chapitres rédigés du niveau).

import type { Contenu, Niveau } from '../contenu/types'
import { melanger } from '../lib/hasard'

export type GenreQuestion = 'vrai-faux' | 'completer' | 'calcul' | 'question'

export interface QuestionQcm {
  id: string
  genre: GenreQuestion
  /** Chapitre et notion (carte) dont vient la question. */
  chapitre: string
  carte: string
  /** Énoncé de l'exercice, pour une question de calcul. */
  contexte?: string
  enonce: string
  propositions: string[]
  bonne: number
  /** Justification, contre-exemple ou correction, montrés après la réponse. */
  explication?: string
  doute?: string
}

export const CONSIGNES: Record<GenreQuestion, string> = {
  'vrai-faux': 'Vrai ou faux ?',
  completer: 'Compléter l’énoncé',
  calcul: 'Question d’exercice',
  question: 'Question',
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
        chapitre,
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
          chapitre,
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
        chapitre,
        carte: e.carte,
        contexte: e.enonce,
        enonce: e.reponse_finale.question,
        propositions: e.reponse_finale.propositions,
        bonne: e.reponse_finale.bonne,
        explication: e.correction,
      })
  for (const q of (c.qcm ?? []).filter((fq) => fq.chapitre === chapitre).flatMap((fq) => fq.questions))
    sortie.push({
      id: q.id,
      genre: q.genre === 'vrai-faux' ? 'vrai-faux' : 'question',
      chapitre,
      carte: q.carte,
      enonce: q.enonce,
      propositions: q.genre === 'vrai-faux' ? ['Vrai', 'Faux'] : q.propositions!,
      bonne: q.genre === 'vrai-faux' ? (q.vrai ? 0 : 1) : q.bonne!,
      explication: q.explication,
      doute: q.doute,
    })
  return sortie
}

/** Mélange les propositions d'une question à choix (« Vrai / Faux » garde son ordre) en suivant la bonne réponse. */
function melangerPropositions(q: QuestionQcm, alea: () => number): QuestionQcm {
  if (q.genre === 'vrai-faux') return q
  const ordre = melanger(q.propositions.map((_, i) => i), alea)
  return { ...q, propositions: ordre.map((i) => q.propositions[i]), bonne: ordre.indexOf(q.bonne) }
}

/** Prend une question dans chaque groupe à tour de rôle (groupes et questions mélangés) jusqu'à en avoir `n`. */
function tourDeRole(groupes: QuestionQcm[][], n: number, alea: () => number): QuestionQcm[] {
  const listes = melanger(groupes.map((g) => melanger(g, alea)), alea)
  const choisies: QuestionQcm[] = []
  for (let tour = 0; choisies.length < n && listes.some((l) => l.length > tour); tour++)
    for (const l of listes) if (l[tour] && choisies.length < n) choisies.push(l[tour])
  return choisies
}

/**
 * Tire `n` questions d'un chapitre en alternant les genres (pour ne pas avoir que des vrai-faux), dans un ordre mélangé,
 * propositions mélangées sauf « Vrai / Faux ».
 */
export function tirerQcm(questions: QuestionQcm[], n = 10, alea: () => number = Math.random): QuestionQcm[] {
  const genres: GenreQuestion[] = ['vrai-faux', 'completer', 'calcul', 'question']
  const choisies = tourDeRole(genres.map((g) => questions.filter((q) => q.genre === g)), n, alea)
  return melanger(choisies, alea).map((q) => melangerPropositions(q, alea))
}

/** Ce que vise un QCM : un chapitre, ou tout un niveau en questions à choix (« qcm ») ou en vrai-faux. */
export type CibleQcm = { chapitre: string } | { niveau: Niveau; sorte: 'qcm' | 'vrai-faux' }

/** Clé sous laquelle les scores sont gardés (progression.qcm) : l'identifiant du chapitre, ou `niveau-<niveau>-<sorte>`. */
export const cleQcm = (c: CibleQcm): string => ('chapitre' in c ? c.chapitre : `niveau-${c.niveau}-${c.sorte}`)

export const NB_QUESTIONS = { chapitre: 10, niveau: 20 }

/** Questions d'un niveau entier : vrai-faux seuls, ou toutes les questions à choix (énoncés à trous, calculs, questions). */
export function questionsNiveau(c: Contenu, niveau: Niveau, sorte: 'qcm' | 'vrai-faux'): QuestionQcm[] {
  return c.themes
    .flatMap((t) => t.chapitres)
    .filter((ch) => ch.niveau === niveau)
    .flatMap((ch) => questionsChapitre(c, ch.id))
    .filter((q) => (sorte === 'vrai-faux') === (q.genre === 'vrai-faux'))
}

/** Questions possibles pour une cible. */
export const questionsCible = (c: Contenu, cible: CibleQcm): QuestionQcm[] =>
  'chapitre' in cible ? questionsChapitre(c, cible.chapitre) : questionsNiveau(c, cible.niveau, cible.sorte)

/**
 * Tirage pour une cible : 10 questions de genres alternés pour un chapitre ; 20 pour un niveau, prises à tour de rôle dans
 * chacun de ses chapitres pour couvrir tout le niveau.
 */
export function tirerCible(c: Contenu, cible: CibleQcm, alea: () => number = Math.random): QuestionQcm[] {
  const questions = questionsCible(c, cible)
  if ('chapitre' in cible) return tirerQcm(questions, NB_QUESTIONS.chapitre, alea)
  const parChapitre = [...new Set(questions.map((q) => q.chapitre))].map((ch) => questions.filter((q) => q.chapitre === ch))
  return melanger(tourDeRole(parChapitre, NB_QUESTIONS.niveau, alea), alea).map((q) => melangerPropositions(q, alea))
}
