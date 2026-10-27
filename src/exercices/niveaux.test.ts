import { describe, expect, it } from 'vitest'
import type { Tentative } from '../lib/progression'
import { exercice } from '../repetition/outils-test'
import { etatExercice, exerciceSuggere, niveauxDe } from './niveaux'

const t = (resultat: Tentative['resultat']): Tentative => ({ date: 0, jour: '2026-03-01', resultat, indices: 0 })
const EX = [exercice('a-ex1', 'a', 1), exercice('a-ex2', 'a', 2), exercice('a-ex3', 'a', 2), exercice('a-ex4', 'a', 4)]

describe('niveaux d’exercices', () => {
  it('grise les niveaux sans exercice (notion sans exercice d’un niveau donné)', () => {
    expect(niveauxDe(EX, {}).map((n) => n.etat)).toEqual(['a-faire', 'a-faire', 'aucun', 'a-faire', 'aucun'])
  })
  it('réussi si tous les exercices du niveau ont été réussis une fois, en partie dès une tentative', () => {
    const h = { 'a-ex1': [t('rate'), t('reussi')], 'a-ex2': [t('reussi')], 'a-ex4': [t('partiel')] }
    expect(niveauxDe(EX, h).map((n) => n.etat)).toEqual(['reussi', 'partiel', 'aucun', 'partiel', 'aucun'])
    expect(etatExercice(h['a-ex1'])).toBe('reussi')
    expect(etatExercice([t('rate'), t('partiel')])).toBe('partiel')
  })
  it('suggère le premier niveau non réussi, puis son premier exercice non réussi', () => {
    expect(exerciceSuggere(EX, {})?.id).toBe('a-ex1')
    expect(exerciceSuggere(EX, { 'a-ex1': [t('reussi')] })?.id).toBe('a-ex2')
    expect(exerciceSuggere(EX, { 'a-ex1': [t('reussi')], 'a-ex2': [t('reussi')] })?.id).toBe('a-ex3')
    const tout = Object.fromEntries(EX.map((e) => [e.id, [t('reussi')]]))
    expect(exerciceSuggere(EX, tout)).toBeNull()
    expect(exerciceSuggere([], {})).toBeNull()
  })
})

describe('variantes', () => {
  it('propose d’abord une variante jamais faite, puis la plus anciennement faite', async () => {
    const { varianteSuivante } = await import('./niveaux')
    expect(varianteSuivante(1, undefined)).toBe(0)
    expect(varianteSuivante(3, [])).toBe(0)
    expect(varianteSuivante(3, [{ ...t('rate'), variante: 0 }])).toBe(1)
    expect(varianteSuivante(3, [t('rate'), { ...t('reussi'), variante: 1 }])).toBe(2) // sans numéro : variante 0
    expect(varianteSuivante(3, [{ ...t('rate'), variante: 2 }, { ...t('rate'), variante: 0 }, { ...t('rate'), variante: 1 }])).toBe(2)
    expect(varianteSuivante(2, [{ ...t('rate'), variante: 5 }])).toBe(0) // variante disparue du contenu : ignorée
  })
})
