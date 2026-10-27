import { describe, expect, it } from 'vitest'
import { avancer, nouvelleSeance, position } from './deroulement'
import type { ElementSeance } from './seance'

const T0 = Date.parse('2026-03-01T09:00:00Z')
const MIN = 60_000
const el = (carte: string): ElementSeance => ({ carte, volets: ['enonce'], neuve: false })

describe('déroulé d’une séance courte', () => {
  it('une carte « À revoir » revient en fin de séance, au plus deux fois', () => {
    let s = nouvelleSeance([el('a'), el('b')], T0)
    s = avancer(s, 'enonce', 0, T0 + MIN, 15)
    expect(s.elements.map((e) => e.carte)).toEqual(['a', 'b', 'a'])
    expect(position(s)).toEqual({ n: 2, total: 2 })
    s = avancer(s, 'enonce', 2, T0 + 2 * MIN, 15) // b
    expect(position(s)).toEqual({ n: 2, total: 2 }) // reprise de a : pas une carte de plus
    s = avancer(s, 'enonce', 0, T0 + 3 * MIN, 15) // a, 2e passage
    s = avancer(s, 'enonce', 0, T0 + 4 * MIN, 15) // a, 3e passage : plus de reprise
    expect(s.fin).toBe(true)
    expect(s.elements.filter((e) => e.carte === 'a')).toHaveLength(3)
    expect(s.vues).toEqual(['a', 'b'])
    expect(s.notes).toEqual([0, 2, 0, 0])
  })
  it('s’arrête quand la durée est écoulée, après la carte en cours', () => {
    let s = nouvelleSeance([el('a'), el('b'), el('c')], T0)
    s = avancer(s, 'enonce', 2, T0 + 5 * MIN, 15)
    expect(s.fin).toBe(false)
    s = avancer(s, 'enonce', 2, T0 + 15 * MIN, 15)
    expect(s.fin).toBe(true)
    expect(s.index).toBe(2)
  })
  it('une séance vide est finie d’emblée, et une séance finie n’avance plus', () => {
    const s = nouvelleSeance([], T0)
    expect(s.fin).toBe(true)
    expect(avancer(s, 'enonce', 2, T0, 15)).toBe(s)
  })
})
