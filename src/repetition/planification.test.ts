import { describe, expect, it } from 'vitest'
import { ajouterJours, jourLocal } from '../lib/temps'
import { delaiAffiche, estDu, planifier, SEUIL_ACQUIS, type EtatVolet, type Note } from './planification'

// Faux temps : 1er mars 2026, 9 h UTC ; jours calculés en UTC.
const T0 = Date.parse('2026-03-01T09:00:00Z')
const H = 3_600_000
const jour = (t: number) => jourLocal(t, 0)

/** Révise un volet jour après jour, toujours à son échéance, avec la même note. */
function enchainer(note: Note, fois: number): EtatVolet {
  let e: EtatVolet | undefined
  let t = T0
  for (let i = 0; i < fois; i++) {
    e = planifier(e, note, t, jour(t))
    t = Date.parse(`${e.echeance}T09:00:00Z`)
  }
  return e!
}

describe('planification SM-2 simplifiée', () => {
  it('affiche les délais d’une carte neuve : < 10 min, 1 j, 4 j, 9 j', () => {
    const d = ([0, 1, 2, 3] as Note[]).map((n) => delaiAffiche(undefined, n, T0, jour(T0)))
    expect(d).toEqual(['< 10 min', '1 j', '4 j', '9 j'])
  })

  it('« À revoir » reprend dix minutes plus tard, le même jour', () => {
    const e = planifier(undefined, 0, T0, jour(T0))
    expect(e.intervalle).toBe(0)
    expect(estDu(e, T0 + 5 * 60_000, jour(T0))).toBe(false)
    expect(estDu(e, T0 + 10 * 60_000, jour(T0))).toBe(true)
    expect(e.echecs).toBe(0) // découvrir une carte sans la savoir n'est pas un oubli
  })

  it('les délais croissent strictement de Difficile à Facile, et avec les réussites', () => {
    const e = enchainer(2, 3)
    const delais = ([1, 2, 3] as Note[]).map((n) => planifier(e, n, T0, e.echeance).intervalle)
    expect(delais[0]).toBeLessThan(delais[1])
    expect(delais[1]).toBeLessThan(delais[2])
    expect(enchainer(2, 4).intervalle).toBeGreaterThan(enchainer(2, 3).intervalle)
  })

  it('un volet devient acquis après quelques « Bien »', () => {
    expect(enchainer(2, 3).intervalle).toBeGreaterThanOrEqual(SEUIL_ACQUIS)
    expect(enchainer(1, 3).intervalle).toBeLessThan(SEUIL_ACQUIS)
  })

  it('la facilité reste bornée', () => {
    expect(enchainer(1, 30).facilite).toBeCloseTo(1.3)
    expect(enchainer(3, 30).facilite).toBeCloseTo(3)
    expect(enchainer(3, 30).intervalle).toBe(365)
  })

  it('oublier un volet acquis le remet en apprentissage et compte un échec', () => {
    const e = enchainer(2, 4)
    const oubli = planifier(e, 0, T0, e.echeance)
    expect(oubli.intervalle).toBe(0)
    expect(oubli.echecs).toBe(1)
    expect(oubli.facilite).toBeLessThan(e.facilite)
    // Reprise prudente : « Bien » donne 2 jours, pas 4.
    expect(planifier(oubli, 2, T0, e.echeance).intervalle).toBe(2)
  })

  it('une carte en retard de plusieurs semaines reste due et gagne un intervalle plus long', () => {
    const e = enchainer(2, 2)
    const enRetard = ajouterJours(e.echeance, 35)
    expect(estDu(e, T0, enRetard)).toBe(true)
    const aLHeure = planifier(e, 2, T0, e.echeance).intervalle
    const tard = planifier(e, 2, T0, enRetard).intervalle
    expect(tard).toBeGreaterThan(aLHeure)
    // Mais « Difficile » ne profite pas du retard.
    expect(planifier(e, 1, T0, enRetard).intervalle).toBe(planifier(e, 1, T0, e.echeance).intervalle)
    // L'échéance part du jour de la révision, pas de l'ancienne échéance.
    expect(planifier(e, 2, T0, enRetard).echeance).toBe(ajouterJours(enRetard, tard))
  })

  it('l’échéance est un jour : un changement de fuseau ne la décale pas d’un jour entier', () => {
    const e = planifier(undefined, 1, T0, jourLocal(T0, 60)) // à Paris : 1 j
    expect(e.echeance).toBe('2026-03-02')
    const lendemainNY = Date.parse('2026-03-02T06:00:00Z') // 1 h du matin à New York, le 2 mars
    expect(estDu(e, lendemainNY, jourLocal(lendemainNY, -300))).toBe(true)
    const veilleNY = Date.parse('2026-03-02T04:00:00Z') // 23 h à New York, encore le 1er mars
    expect(estDu(e, veilleNY, jourLocal(veilleNY, -300))).toBe(false)
  })

  it('une reprise « À revoir » de la veille est due le lendemain même avant l’heure', () => {
    const tard = Date.parse('2026-03-01T23:55:00Z')
    const e = planifier(undefined, 0, tard, jour(tard))
    const minuitUne = tard + 6 * 60_000 // le 2 mars à 0 h 01 : nouveau jour, reprise pas encore échue en minutes
    expect(estDu(e, minuitUne, jour(minuitUne))).toBe(true)
  })

  it('compte révisions et réussites', () => {
    let e = planifier(undefined, 0, T0, jour(T0))
    e = planifier(e, 2, T0 + H, jour(T0 + H))
    e = planifier(e, 1, T0 + 2 * H, jour(T0 + 2 * H))
    expect([e.revisions, e.reussites]).toEqual([3, 1])
    expect(e.premiere).toBe('2026-03-01')
  })
})
