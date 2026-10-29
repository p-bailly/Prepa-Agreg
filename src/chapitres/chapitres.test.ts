import katex from 'katex'
import { describe, expect, it } from 'vitest'
import { CONTENU } from '../contenu/chargement'
import { formulesDe } from '../contenu/texte'
import { enregistrerQcm, enregistrerTentative, PROGRESSION_VIDE, apprendre } from '../lib/progression'
import { carte, contenu, exercice } from '../repetition/outils-test'
import { avancementChapitre, avancementNotion, chapitresNumerotes, meilleurQcm } from './plan'
import { cleQcm, questionsChapitre, questionsNiveau, tirerCible, tirerQcm } from './qcm'

const J = '2026-03-10'
const T = Date.parse(`${J}T09:00:00Z`)

/** Tirage pseudo-aléatoire reproductible. */
function graine(s: number) {
  return () => ((s = (s * 16807) % 2147483647) - 1) / 2147483646
}

describe('plan des chapitres', () => {
  it('numérote tous les chapitres de 1 à n, Terminale d’abord', () => {
    const l = chapitresNumerotes(CONTENU)
    expect(l.map((x) => x.numero)).toEqual(l.map((_, i) => i + 1))
    expect(l.length).toBe(CONTENU.themes.reduce((s, t) => s + t.chapitres.length, 0))
    expect(l[0].chapitre.niveau).toBe('terminale')
    expect(l[l.length - 1].chapitre.niveau).toBe('agreg')
  })

  it('suit l’avancement d’une notion et d’un chapitre', () => {
    const c = contenu([carte('a', 'x', 'l1'), carte('b', 'x', 'l1')], [exercice('a-ex1', 'a', 1), exercice('b-ex1', 'b', 2)])
    const exA = c.exercices.filter((e) => e.carte === 'a')
    let p = PROGRESSION_VIDE
    expect(avancementNotion(c.cartes[0], exA, p)).toBe('a-faire')
    expect(avancementChapitre(c, 'x-l1', p)).toBe('a-faire')
    expect(avancementChapitre(c, 'inconnu', p)).toBe('vide')
    p = apprendre(p, ['a'], J)
    expect(avancementNotion(c.cartes[0], exA, p)).toBe('en-cours')
    p = enregistrerTentative(p, 'a-ex1', { date: T, jour: J, resultat: 'reussi', indices: 0 })
    p = enregistrerTentative(p, 'b-ex1', { date: T, jour: J, resultat: 'reussi', indices: 0 })
    expect(avancementNotion(c.cartes[0], exA, p)).toBe('fait')
    expect(avancementChapitre(c, 'x-l1', p)).toBe('en-cours')
    p = enregistrerQcm(p, 'x-l1', { date: T, jour: J, bonnes: 7, total: 10 })
    expect(avancementChapitre(c, 'x-l1', p)).toBe('en-cours')
    p = enregistrerQcm(p, 'x-l1', { date: T, jour: J, bonnes: 9, total: 10 })
    expect(meilleurQcm(p, 'x-l1')).toMatchObject({ bonnes: 9 })
    expect(avancementChapitre(c, 'x-l1', p)).toBe('fait')
  })
})

describe('QCM des chapitres', () => {
  const rediges = [...new Set(CONTENU.cartes.map((k) => k.chapitre))]

  it('offre au moins 10 questions bien formées par chapitre rédigé', () => {
    for (const ch of rediges) {
      const qs = questionsChapitre(CONTENU, ch)
      expect(qs.length, ch).toBeGreaterThanOrEqual(10)
      for (const q of qs) {
        expect(q.enonce, q.id).not.toMatch(/\{\{\d+\}\}/)
        expect(q.bonne, q.id).toBeGreaterThanOrEqual(0)
        expect(q.bonne, q.id).toBeLessThan(q.propositions.length)
        expect(new Set(q.propositions).size, q.id).toBe(q.propositions.length)
        for (const t of [q.enonce, q.contexte, q.explication, ...q.propositions])
          for (const f of t ? formulesDe(t) : []) expect(() => katex.renderToString(f.t, { throwOnError: true, displayMode: f.bloc }), q.id).not.toThrow()
      }
    }
  })

  it('reprend les questions écrites du chapitre', () => {
    const a = carte('a', 'x', 'l1')
    const c = {
      ...contenu([a]),
      qcm: [
        {
          chapitre: 'x-l1',
          questions: [
            { id: 'x-l1-q1', carte: 'a', genre: 'vrai-faux' as const, enonce: 'e', vrai: false, explication: 'x' },
            { id: 'x-l1-q2', carte: 'a', genre: 'choix' as const, enonce: 'e', propositions: ['1', '2', '3', '4'], bonne: 2, explication: 'x' },
          ],
        },
      ],
    }
    const qs = questionsChapitre(c, 'x-l1')
    expect(qs.map((q) => [q.genre, q.propositions[q.bonne]])).toEqual([
      ['vrai-faux', 'Faux'],
      ['question', '3'],
    ])
  })

  it('tire des questions de chaque genre et garde la bonne réponse en mélangeant', () => {
    const qs = questionsChapitre(CONTENU, 'ar-t-suites')
    const t = tirerQcm(qs, 10, graine(7))
    expect(t).toHaveLength(10)
    expect(new Set(t.map((q) => q.id)).size).toBe(10)
    expect(new Set(t.map((q) => q.genre)).size).toBeGreaterThanOrEqual(3)
    for (const q of t) {
      const origine = qs.find((x) => x.id === q.id)!
      expect(q.propositions[q.bonne]).toBe(origine.propositions[origine.bonne])
    }
    expect(tirerQcm(qs, 10, graine(7)).map((q) => q.id)).toEqual(t.map((q) => q.id))
  })

  it('tire un QCM et un vrai-faux de niveau qui couvrent ses chapitres', () => {
    for (const niveau of ['terminale', 'l1'] as const) {
      const chapitres = new Set(CONTENU.cartes.filter((k) => k.niveau === niveau).map((k) => k.chapitre))
      const vf = tirerCible(CONTENU, { niveau, sorte: 'vrai-faux' }, graine(3))
      const qcm = tirerCible(CONTENU, { niveau, sorte: 'qcm' }, graine(5))
      expect(vf).toHaveLength(20)
      expect(qcm).toHaveLength(20)
      expect(vf.every((q) => q.genre === 'vrai-faux')).toBe(true)
      expect(qcm.every((q) => q.genre !== 'vrai-faux' && q.propositions.length === 4)).toBe(true)
      // À tour de rôle : 20 questions pour au plus 20 chapitres, donc min(20, nombre de chapitres) chapitres distincts.
      for (const t of [vf, qcm]) {
        expect(new Set(t.map((q) => q.id)).size).toBe(20)
        expect(new Set(t.map((q) => q.chapitre)).size).toBe(Math.min(20, chapitres.size))
        for (const q of t) expect(chapitres.has(q.chapitre)).toBe(true)
      }
    }
    expect(questionsNiveau(CONTENU, 'agreg', 'qcm')).toEqual([])
    expect(cleQcm({ niveau: 'l1', sorte: 'qcm' })).toBe('niveau-l1-qcm')
    expect(cleQcm({ chapitre: 'ar-t-suites' })).toBe('ar-t-suites')
  })
})
