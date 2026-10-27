import { describe, expect, it } from 'vitest'
import { enregistrerTentative, noterVolet, PROGRESSION_VIDE, type Progression } from '../lib/progression'
import { ajouterJours } from '../lib/temps'
import { carte, contenu, exercice, reglages } from '../repetition/outils-test'
import {
  cartesAcquises,
  cartesParJour,
  exercicesParDifficulte,
  joursAvantEcrit,
  niveauAtteint,
  progresNiveaux,
  reussiteParTheme,
  serieJours,
} from './tableau'

const J = '2026-03-10'
const T = Date.parse(`${J}T09:00:00Z`)

/** Fait passer une carte à « acquise » (trois « Bien » espacés). */
function acquerir(p: Progression, id: string): Progression {
  let q = p
  let j = '2026-01-01'
  for (let k = 0; k < 4; k++) {
    q = noterVolet(q, id, 'enonce', 2, T, j, true)
    j = q.volets[`${id}#enonce`].echeance
  }
  return q
}

describe('tableau de bord', () => {
  it('série de jours : consécutifs, jusqu’à aujourd’hui ou hier', () => {
    let p = PROGRESSION_VIDE
    for (const d of [-3, -2, -1]) p = noterVolet(p, 'x', 'enonce', 2, T, ajouterJours(J, d), true)
    expect(serieJours(p, J)).toBe(3) // rien encore aujourd'hui : la série d'hier tient
    p = noterVolet(p, 'x', 'enonce', 2, T, J, true)
    expect(serieJours(p, J)).toBe(4)
    expect(serieJours(p, ajouterJours(J, 2))).toBe(0) // un jour manqué
    expect(serieJours(PROGRESSION_VIDE, J)).toBe(0)
  })

  it('progression par niveau, rapportée aux cartes visées, et niveau atteint à 80 %', () => {
    const C = contenu([carte('t1', 'a', 'terminale'), carte('t2', 'a', 'terminale'), carte('l1', 'a', 'l1')])
    // contenu() vise exactement les cartes existantes ; on ajoute 3 cartes visées en L1.
    C.themes[0].chapitres[1].prevues = 4
    let p = acquerir(PROGRESSION_VIDE, 't1')
    let n = progresNiveaux(C, p, reglages())
    expect(n.map((x) => [x.acquises, x.visees])).toEqual([[1, 2], [0, 4], [0, 0], [0, 0], [0, 0]])
    expect(niveauAtteint(n)).toBeNull()
    p = acquerir(p, 't2')
    p = acquerir(p, 'l1')
    n = progresNiveaux(C, p, reglages())
    expect(niveauAtteint(n)).toBe('terminale') // L1 : 1 / 4 seulement
    expect(cartesAcquises(C, p, reglages())).toBe(3)
  })

  it('niveau de départ relevé : les niveaux inférieurs comptent comme acquis, ou sont masqués', () => {
    const C = contenu([carte('t1', 'a', 'terminale'), carte('l1', 'a', 'l1')])
    const acquis = progresNiveaux(C, PROGRESSION_VIDE, reglages({ niveauDepart: 'l1' }))
    expect(acquis[0].statut).toBe('acquis')
    expect(niveauAtteint(acquis)).toBe('terminale')
    expect(cartesAcquises(C, PROGRESSION_VIDE, reglages({ niveauDepart: 'l1' }))).toBe(1)
    const masques = reglages({ niveauDepart: 'l1', niveauxInferieurs: 'masques' })
    expect(niveauAtteint(progresNiveaux(C, PROGRESSION_VIDE, masques))).toBeNull()
    expect(cartesAcquises(C, PROGRESSION_VIDE, masques)).toBe(0)
  })

  it('histogramme des 30 derniers jours, aujourd’hui en dernier', () => {
    let p = noterVolet(PROGRESSION_VIDE, 'x', 'enonce', 2, T, J, true)
    p = noterVolet(p, 'y', 'enonce', 2, T, J, true)
    p = noterVolet(p, 'x', 'demo', 2, T, J, false) // un volet secondaire ne compte pas une carte de plus
    p = noterVolet(p, 'x', 'enonce', 2, T, ajouterJours(J, -29), true)
    p = noterVolet(p, 'x', 'enonce', 2, T, ajouterJours(J, -30), true) // hors fenêtre
    const h = cartesParJour(p, J)
    expect(h).toHaveLength(30)
    expect(h[29]).toEqual({ jour: J, cartes: 2 })
    expect(h[0].cartes).toBe(1)
    expect(h.reduce((s, x) => s + x.cartes, 0)).toBe(3)
  })

  it('réussite par thème, faibles en premier', () => {
    const C = contenu([carte('a1', 'analyse', 'l1'), carte('b1', 'algebre', 'l1')])
    let p = noterVolet(PROGRESSION_VIDE, 'a1', 'enonce', 2, T, J, true)
    p = noterVolet(p, 'b1', 'enonce', 0, T, J, true)
    p = noterVolet(p, 'b1', 'enonce', 3, T, J, true)
    expect(reussiteParTheme(C, p).map((x) => [x.theme, x.part])).toEqual([
      ['algebre', 0.5],
      ['analyse', 1],
    ])
  })

  it('exercices réussis par difficulté', () => {
    const C = contenu([carte('a', 't', 'l1')], [exercice('a-ex1', 'a', 1), exercice('a-ex2', 'a', 3), exercice('a-ex3', 'a', 3)])
    let p = enregistrerTentative(PROGRESSION_VIDE, 'a-ex2', { date: T, jour: J, resultat: 'reussi', indices: 1 })
    p = enregistrerTentative(p, 'a-ex3', { date: T, jour: J, resultat: 'partiel', indices: 2 })
    expect(exercicesParDifficulte(C, p)).toEqual([
      { difficulte: 1, reussis: 0, total: 1 },
      { difficulte: 2, reussis: 0, total: 0 },
      { difficulte: 3, reussis: 1, total: 2 },
      { difficulte: 4, reussis: 0, total: 0 },
      { difficulte: 5, reussis: 0, total: 0 },
    ])
    expect(p.jours[J].exercices).toBe(2)
  })

  it('compte à rebours seulement avec une date', () => {
    expect(joursAvantEcrit(null, J)).toBeNull()
    expect(joursAvantEcrit('2027-03-04', J)).toBe(359)
  })
})
