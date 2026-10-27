import { describe, expect, it, beforeEach } from 'vitest'
import { ecrire, lire } from './stockage'

// localStorage minimal pour Node.
const memoire = new Map<string, string>()
globalThis.localStorage = {
  getItem: (c: string) => memoire.get(c) ?? null,
  setItem: (c: string, v: string) => void memoire.set(c, v),
  removeItem: (c: string) => void memoire.delete(c),
  clear: () => memoire.clear(),
  key: () => null,
  length: 0,
} as Storage

describe('stockage', () => {
  beforeEach(() => memoire.clear())
  it('fusionne la valeur lue avec les valeurs par défaut', () => {
    ecrire('reglages-v1', { a: 2 })
    expect(lire('reglages-v1', { a: 1, b: 3 })).toEqual({ a: 2, b: 3 })
  })
  it('renvoie la valeur par défaut si la donnée est illisible', () => {
    memoire.set('x', '{pas du json')
    expect(lire('x', { a: 1 })).toEqual({ a: 1 })
    memoire.set('x', '[1,2]')
    expect(lire('x', { a: 1 })).toEqual({ a: 1 })
  })
})
