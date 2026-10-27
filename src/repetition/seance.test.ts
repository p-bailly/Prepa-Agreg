import { describe, expect, it } from 'vitest'
import { PROGRESSION_VIDE, noterVolet, type Progression } from '../lib/progression'
import { ajouterJours, jourLocal } from '../lib/temps'
import { etatAffiche, etatCarte, prerequisManquants, voletsDus } from './etat'
import { choisirMode, modesPossibles } from './modes'
import { carte, contenu, reglages } from './outils-test'
import { cartesARevoir, cartesNeuves, construireSeance, neuvesDuJour, ordreIntroduction, resumeDuJour } from './seance'

const T0 = Date.parse('2026-03-01T09:00:00Z')
const J0 = jourLocal(T0, 0)

/** Note « Bien » tous les volets présents d'une carte. */
function voir(p: Progression, id: string, volets: ('enonce' | 'demo' | 'contre')[] = ['enonce'], t = T0, j = J0) {
  return volets.reduce((q, v, i) => noterVolet(q, id, v, 2, t, j, i === 0), p)
}

describe('ordre d’introduction et prérequis', () => {
  const C = contenu([
    carte('a1', 'analyse', 'l1'),
    carte('a0', 'analyse', 'terminale'),
    carte('b0', 'algebre', 'terminale'),
    carte('a0bis', 'analyse', 'terminale', ['a0']),
  ])
  it('niveau le plus bas d’abord, thèmes alternés, ordre du fichier', () => {
    expect(ordreIntroduction(C).map((c) => c.id)).toEqual(['a0', 'b0', 'a0bis', 'a1'])
  })
  it('une carte neuve attend que ses prérequis soient au moins en cours', () => {
    const parId = new Map(C.cartes.map((c) => [c.id, c]))
    const c = parId.get('a0bis')!
    expect(prerequisManquants(c, parId, {}, reglages())).toEqual(['a0'])
    const p = voir(PROGRESSION_VIDE, 'a0')
    expect(prerequisManquants(c, parId, p.volets, reglages())).toEqual([])
  })
  it('dans une même séance, une carte introduite débloque celles qui en dépendent', () => {
    expect(cartesNeuves(C, PROGRESSION_VIDE, reglages(), 10).map((e) => e.carte)).toEqual(['a0', 'b0', 'a0bis', 'a1'])
  })
  it('un prérequis jamais vu bloque, même s’il est plus loin dans l’ordre', () => {
    const D = contenu([carte('x', 'analyse', 'terminale', ['y']), carte('y', 'analyse', 'terminale')])
    // y arrive après x dans le fichier : x entre à la passe suivante, après y.
    expect(cartesNeuves(D, PROGRESSION_VIDE, reglages(), 10).map((e) => e.carte)).toEqual(['y', 'x'])
    expect(cartesNeuves(D, PROGRESSION_VIDE, reglages(), 1).map((e) => e.carte)).toEqual(['y'])
    // Prérequis inexistant : la carte n'est jamais proposée.
    const E = contenu([carte('z', 'analyse', 'terminale', ['fantome'])])
    expect(cartesNeuves(E, PROGRESSION_VIDE, reglages(), 10)).toEqual([])
  })
  it('un niveau de départ plus haut lève les prérequis des niveaux inférieurs et retire leurs cartes', () => {
    const D = contenu([carte('t', 'analyse', 'terminale'), carte('l', 'analyse', 'l1', ['t'])])
    const r = reglages({ niveauDepart: 'l1' })
    expect(cartesNeuves(D, PROGRESSION_VIDE, r, 10).map((e) => e.carte)).toEqual(['l'])
    expect(etatAffiche(D.cartes[0], {}, r)).toBe('acquise')
    expect(etatAffiche(D.cartes[0], {}, { ...r, niveauxInferieurs: 'masques' })).toBe('neuve')
    // Revenir au niveau de départ Terminale : la carte de Terminale redevient à introduire, en premier.
    expect(cartesNeuves(D, PROGRESSION_VIDE, reglages(), 10).map((e) => e.carte)).toEqual(['t', 'l'])
  })
})

describe('plafond de cartes neuves et séance', () => {
  const C = contenu(Array.from({ length: 30 }, (_, i) => carte(`c${String(i).padStart(2, '0')}`, 'analyse', 'terminale')))
  it('respecte le plafond quotidien, même sur plusieurs séances le même jour', () => {
    const r = reglages({ plafondNeuves: 10, seanceCartes: 6 })
    const s1 = construireSeance(C, PROGRESSION_VIDE, r, T0, J0)
    expect(s1).toHaveLength(6)
    let p = s1.reduce((q, e) => voir(q, e.carte), PROGRESSION_VIDE)
    expect(neuvesDuJour(C, p, J0)).toBe(6)
    const s2 = construireSeance(C, p, r, T0, J0)
    expect(s2.filter((e) => e.neuve)).toHaveLength(4) // plafond atteint : 6 + 4 = 10
    p = s2.reduce((q, e) => voir(q, e.carte), p)
    expect(construireSeance(C, p, r, T0, J0)).toEqual([])
    // Le lendemain, le plafond repart à zéro.
    const T1 = T0 + 86_400_000
    expect(cartesNeuves(C, p, r, 100)).toHaveLength(20)
    expect(construireSeance(C, p, r, T1, ajouterJours(J0, 1)).filter((e) => e.neuve).length).toBeGreaterThan(0)
  })
  it('plafond à zéro : aucune carte neuve', () => {
    expect(construireSeance(C, PROGRESSION_VIDE, reglages({ plafondNeuves: 0 }), T0, J0)).toEqual([])
  })
  it('les révisions passent d’abord, les neuves sont réparties entre elles', () => {
    let p = PROGRESSION_VIDE
    for (const id of ['c00', 'c01', 'c02', 'c03']) p = noterVolet(p, id, 'enonce', 1, T0, J0, true) // 1 j
    const J1 = ajouterJours(J0, 1)
    const s = construireSeance(C, p, reglages({ seanceCartes: 8, plafondNeuves: 4 }), T0 + 86_400_000, J1)
    expect(s.filter((e) => !e.neuve).map((e) => e.carte)).toEqual(['c00', 'c01', 'c02', 'c03'])
    expect(s.filter((e) => e.neuve)).toHaveLength(4)
    expect(s[0].neuve).toBe(false)
    expect(s.map((e) => e.neuve)).toEqual([false, true, false, true, false, true, false, true])
  })
  it('les cartes les plus en retard sortent en premier, même après plusieurs semaines', () => {
    let p = noterVolet(PROGRESSION_VIDE, 'c00', 'enonce', 1, T0, J0, true) // due J0+1
    p = noterVolet(p, 'c01', 'enonce', 3, T0, J0, true) // due J0+9
    const plusTard = ajouterJours(J0, 40)
    expect(cartesARevoir(C, p, reglages(), T0, plusTard).map((e) => e.carte)).toEqual(['c00', 'c01'])
  })
  it('résumé du jour : à revoir, nouvelles, minutes estimées bornées par la durée de séance', () => {
    const res = resumeDuJour(C, PROGRESSION_VIDE, reglages(), T0, J0)
    expect(res).toEqual({ aRevoir: 0, nouvelles: 10, minutes: 15 })
  })
  it('filtres par thème, niveau et état', () => {
    const D = contenu([carte('a', 'analyse', 'terminale'), carte('b', 'algebre', 'terminale'), carte('c', 'algebre', 'l1')])
    const ids = (r: ReturnType<typeof reglages>) => construireSeance(D, PROGRESSION_VIDE, r, T0, J0).map((e) => e.carte)
    expect(ids(reglages({ filtres: { themes: ['algebre'], niveaux: [], etats: [] } }))).toEqual(['b', 'c'])
    expect(ids(reglages({ filtres: { themes: [], niveaux: ['l1'], etats: [] } }))).toEqual(['c'])
    expect(ids(reglages({ filtres: { themes: [], niveaux: [], etats: ['en-cours'] } }))).toEqual([])
    expect(ids(reglages({ niveauxActifs: ['terminale'] }))).toEqual(['a', 'b'])
  })
})

describe('volets', () => {
  const c = carte('t', 'analyse', 'l1', [], { demo: 'Preuve.', contre_exemple: 'Contre-exemple.' })
  const C = contenu([c])
  it('chaque volet a son propre état ; la carte n’est acquise que si tous le sont', () => {
    let p = voir(PROGRESSION_VIDE, 't', ['enonce'])
    expect(etatCarte(c, p.volets)).toBe('en-cours')
    // Démonstration et contre-exemple jamais ouverts : à voir dès la séance suivante.
    expect(voletsDus(c, p.volets, T0, J0)).toEqual(['demo', 'contre'])
    p = voir(p, 't', ['demo', 'contre'])
    expect(voletsDus(c, p.volets, T0, J0)).toEqual([])
    // Trois semaines de « Bien » sur tout : acquise.
    let t = T0
    let j = J0
    for (let k = 0; k < 6 && etatCarte(c, p.volets) !== 'acquise'; k++) {
      const dus = voletsDus(c, p.volets, t, j)
      p = voir(p, 't', dus.length ? dus : [], t, j)
      j = ajouterJours(j, 30)
      t += 30 * 86_400_000
    }
    expect(etatCarte(c, p.volets)).toBe('acquise')
    expect(cartesARevoir(C, p, reglages(), T0, J0)).toEqual([])
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
    // Des affirmations toutes vraies ne permettent pas le vrai/faux avec contre-exemple.
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
