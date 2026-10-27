import { describe, expect, it } from 'vitest'
import { apprendre, noterVolet, PROGRESSION_VIDE, type Progression } from '../lib/progression'
import { ajouterJours, jourLocal } from '../lib/temps'
import { etatAffiche, etatCarte, prerequisManquants, voletsDus } from './etat'
import { choisirMode, modesPossibles } from './modes'
import { carte, contenu, reglages } from './outils-test'
import { aApprendre, cartesARevoir, construireSeance, leconDuJour, neuvesDuJour, ordreIntroduction, resumeDuJour } from './seance'

const T0 = Date.parse('2026-03-01T09:00:00Z')
const J0 = jourLocal(T0, 0)
const J1 = ajouterJours(J0, 1)
const T1 = T0 + 86_400_000
const VIDE = PROGRESSION_VIDE

/** Note « Bien » des volets d'une carte. */
function voir(p: Progression, id: string, volets: ('enonce' | 'demo' | 'contre')[] = ['enonce'], t = T0, j = J0) {
  return volets.reduce((q, v, i) => noterVolet(q, id, v, 2, t, j, i === 0), p)
}
const ids = (l: { id: string }[]) => l.map((c) => c.id)

describe('ordre d’introduction et prérequis', () => {
  const C = contenu([
    carte('a1', 'analyse', 'l1'),
    carte('a0', 'analyse', 'terminale'),
    carte('b0', 'algebre', 'terminale'),
    carte('a0bis', 'analyse', 'terminale', ['a0']),
  ])
  it('niveau le plus bas d’abord, thèmes alternés, ordre du fichier', () => {
    expect(ids(ordreIntroduction(C))).toEqual(['a0', 'b0', 'a0bis', 'a1'])
  })
  it('une notion attend que ses prérequis soient au moins appris', () => {
    const parId = new Map(C.cartes.map((c) => [c.id, c]))
    const c = parId.get('a0bis')!
    expect(prerequisManquants(c, parId, VIDE, reglages())).toEqual(['a0'])
    expect(prerequisManquants(c, parId, apprendre(VIDE, ['a0'], J0), reglages())).toEqual([])
  })
  it('dans une même leçon, une notion apprise débloque celles qui en dépendent', () => {
    expect(ids(aApprendre(C, VIDE, reglages(), 10))).toEqual(['a0', 'b0', 'a0bis', 'a1'])
  })
  it('un prérequis jamais vu bloque, même s’il est plus loin dans l’ordre', () => {
    const D = contenu([carte('x', 'analyse', 'terminale', ['y']), carte('y', 'analyse', 'terminale')])
    expect(ids(aApprendre(D, VIDE, reglages(), 10))).toEqual(['y', 'x'])
    expect(ids(aApprendre(D, VIDE, reglages(), 1))).toEqual(['y'])
    const E = contenu([carte('z', 'analyse', 'terminale', ['fantome'])])
    expect(aApprendre(E, VIDE, reglages(), 10)).toEqual([])
  })
  it('un niveau de départ plus haut lève les prérequis des niveaux inférieurs et retire leurs notions', () => {
    const D = contenu([carte('t', 'analyse', 'terminale'), carte('l', 'analyse', 'l1', ['t'])])
    const r = reglages({ niveauDepart: 'l1' })
    expect(ids(aApprendre(D, VIDE, r, 10))).toEqual(['l'])
    expect(etatAffiche(D.cartes[0], VIDE, r)).toBe('acquise')
    expect(etatAffiche(D.cartes[0], VIDE, { ...r, niveauxInferieurs: 'masques' })).toBe('neuve')
    // Revenir au niveau de départ Terminale : la notion de Terminale redevient à apprendre, en premier.
    expect(ids(aApprendre(D, VIDE, reglages(), 10))).toEqual(['t', 'l'])
  })
})

describe('leçon d’abord, révision ensuite', () => {
  const C = contenu([
    carte('s1', 'analyse', 'terminale', [], { chapitre: 'suites' }),
    carte('s2', 'analyse', 'terminale', [], { chapitre: 'suites' }),
    carte('c1', 'analyse', 'terminale', [], { chapitre: 'continuite' }),
  ])
  C.themes[0].chapitres = [
    { id: 'suites', niveau: 'terminale', titre: 'Suites', notions: ['s1', 's2'], prevues: 2 },
    { id: 'continuite', niveau: 'terminale', titre: 'Continuité', notions: ['c1'], prevues: 1 },
  ]
  it('la leçon du jour regroupe les prochaines notions d’un même chapitre', () => {
    const l = leconDuJour(C, VIDE, reglages(), J0)!
    expect(l.chapitre.id).toBe('suites')
    expect(ids(l.cartes)).toEqual(['s1', 's2'])
    // Le chapitre suivant vient dans une leçon suivante, le même jour si le plafond le permet.
    const apres = apprendre(VIDE, ['s1', 's2'], J0)
    expect(ids(leconDuJour(C, apres, reglages(), J0)!.cartes)).toEqual(['c1'])
    expect(leconDuJour(C, apres, reglages({ plafondNeuves: 2 }), J0)).toBeNull()
  })
  it('jamais de révision d’une notion pas encore apprise, ni le jour même de la leçon', () => {
    expect(construireSeance(C, VIDE, reglages(), T0, J0)).toEqual([])
    const p = apprendre(VIDE, ['s1'], J0)
    expect(etatCarte(C.cartes[0], p)).toBe('en-cours')
    expect(construireSeance(C, p, reglages(), T0, J0)).toEqual([])
    // Le lendemain : première révision de la notion apprise.
    expect(construireSeance(C, p, reglages(), T1, J1).map((e) => e.carte)).toEqual(['s1'])
    expect(resumeDuJour(C, p, reglages(), T1, J1).aRevoir).toBe(1)
  })
  it('le plafond quotidien compte les notions apprises aujourd’hui', () => {
    const r = reglages({ plafondNeuves: 2 })
    expect(resumeDuJour(C, VIDE, r, T0, J0).aApprendre).toBe(2)
    const p = apprendre(VIDE, ['s1'], J0)
    expect(neuvesDuJour(C, p, J0)).toBe(1)
    expect(resumeDuJour(C, p, r, T0, J0).aApprendre).toBe(1)
    // Le lendemain, le plafond repart à zéro.
    expect(resumeDuJour(C, p, r, T1, J1).aApprendre).toBe(2)
    expect(resumeDuJour(C, VIDE, reglages({ plafondNeuves: 0 }), T0, J0).lecon).toBeNull()
  })
  it('apprendre deux fois ne change pas le jour d’apprentissage', () => {
    const p = apprendre(apprendre(VIDE, ['s1'], J0), ['s1', 's2'], J1)
    expect(p.apprises).toEqual({ s1: J0, s2: J1 })
  })
})

describe('séance de révision', () => {
  const C = contenu(Array.from({ length: 30 }, (_, i) => carte(`c${String(i).padStart(2, '0')}`, 'analyse', 'terminale')))
  it('au plus N cartes, les plus en retard d’abord, même après plusieurs semaines', () => {
    let p = noterVolet(VIDE, 'c00', 'enonce', 1, T0, J0, true) // due J0+1
    p = noterVolet(p, 'c01', 'enonce', 3, T0, J0, true) // due J0+9
    const plusTard = ajouterJours(J0, 40)
    expect(cartesARevoir(C, p, reglages(), T0, plusTard).map((e) => e.carte)).toEqual(['c00', 'c01'])
    const tout = apprendre(VIDE, C.cartes.map((c) => c.id), J0)
    expect(construireSeance(C, tout, reglages({ seanceCartes: 6 }), T1, J1)).toHaveLength(6)
  })
  it('résumé du jour : révisions et leçon, minutes bornées par la durée de séance', () => {
    const r = resumeDuJour(C, VIDE, reglages(), T0, J0)
    expect([r.aRevoir, r.aApprendre, r.lecon?.cartes.length, r.minutesRevision]).toEqual([0, 10, 10, 0])
    const tout = apprendre(VIDE, C.cartes.map((c) => c.id), J0)
    expect(resumeDuJour(C, tout, reglages(), T1, J1).minutesRevision).toBe(12) // 15 cartes × 0,75 min
  })
  it('filtres par thème, niveau et état', () => {
    const D = contenu([carte('a', 'analyse', 'terminale'), carte('b', 'algebre', 'terminale'), carte('c', 'algebre', 'l1')])
    const p = apprendre(VIDE, ['a', 'b', 'c'], J0)
    const s = (r: ReturnType<typeof reglages>) => construireSeance(D, p, r, T1, J1).map((e) => e.carte)
    expect(s(reglages({ filtres: { themes: ['algebre'], niveaux: [], etats: [] } }))).toEqual(['b', 'c'])
    expect(s(reglages({ filtres: { themes: [], niveaux: ['l1'], etats: [] } }))).toEqual(['c'])
    expect(s(reglages({ filtres: { themes: [], niveaux: [], etats: ['acquise'] } }))).toEqual([])
    expect(s(reglages({ niveauxActifs: ['terminale'] }))).toEqual(['a', 'b'])
    expect(ids(aApprendre(D, VIDE, reglages({ filtres: { themes: ['algebre'], niveaux: [], etats: [] } }), 10))).toEqual(['b', 'c'])
  })
})

describe('volets', () => {
  const c = carte('t', 'analyse', 'l1', [], { demo: 'Preuve.', contre_exemple: 'Contre-exemple.' })
  const C = contenu([c])
  it('chaque volet a son propre état ; la carte n’est acquise que si tous le sont', () => {
    let p = apprendre(VIDE, ['t'], J0)
    expect(voletsDus(c, p, T0, J0)).toEqual([])
    expect(voletsDus(c, p, T1, J1)).toEqual(['enonce', 'demo', 'contre'])
    p = voir(p, 't', ['enonce'], T1, J1)
    // Démonstration et contre-exemple pas encore notés : toujours dus.
    expect(voletsDus(c, p, T1, J1)).toEqual(['demo', 'contre'])
    p = voir(p, 't', ['demo', 'contre'], T1, J1)
    expect(voletsDus(c, p, T1, J1)).toEqual([])
    let t = T1
    let j = J1
    for (let k = 0; k < 8 && etatCarte(c, p) !== 'acquise'; k++) {
      j = ajouterJours(j, 30)
      t += 30 * 86_400_000
      p = voir(p, 't', voletsDus(c, p, t, j), t, j)
    }
    expect(etatCarte(c, p)).toBe('acquise')
    expect(cartesARevoir(C, p, reglages(), T1, J1)).toEqual([])
  })
})

describe('modes de révision', () => {
  const pleine = carte('p', 'analyse', 'l1', [], {
    hypotheses: ['h'],
    hypotheses_pieges: ['piège'],
    enonce_trous: 'f est {{1}}',
    trous: [{ reponse: 'continue', leurres: ['x', 'y', 'z'] }],
    affirmations: [{ texte: 'faux', vrai: false, contre_exemple: 'ce' }],
  })
  it('un mode n’apparaît que si la carte a ce qu’il faut', () => {
    expect(modesPossibles(carte('n', 'a', 'l1'))).toEqual(['tete'])
    expect(modesPossibles(pleine)).toEqual(['tete', 'hypotheses', 'trous', 'vraifaux'])
    expect(modesPossibles(carte('v', 'a', 'l1', [], { affirmations: [{ texte: 't', vrai: true }] }))).toEqual(['tete'])
  })
  it('alterne d’une carte à l’autre parmi les modes actifs', () => {
    const tous = ['tete', 'hypotheses', 'trous', 'vraifaux'] as const
    expect(choisirMode(pleine, tous, null, false)).toBe('tete')
    expect(choisirMode(pleine, tous, 'tete', false)).toBe('hypotheses')
    expect(choisirMode(pleine, tous, 'vraifaux', false)).toBe('tete')
    expect(choisirMode(pleine, ['trous', 'vraifaux'], 'trous', false)).toBe('vraifaux')
    expect(choisirMode(pleine, ['hypotheses'], 'hypotheses', false)).toBe('hypotheses')
    expect(choisirMode(carte('n', 'a', 'l1'), ['trous'], null, false)).toBe('tete')
    expect(choisirMode(pleine, tous, 'tete', true)).toBe('tete')
  })
})
