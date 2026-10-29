import { describe, expect, it } from 'vitest'
import { ajouterJours, jourLocal, joursAvantEcrit, joursEntre } from './temps'

describe('jours', () => {
  it('donne le jour local selon le fuseau', () => {
    const t = Date.parse('2026-03-01T23:30:00Z')
    expect(jourLocal(t, 0)).toBe('2026-03-01')
    expect(jourLocal(t, 60)).toBe('2026-03-02') // Paris : déjà le lendemain
    expect(jourLocal(t, -300)).toBe('2026-03-01') // New York
  })
  it('compte les jours, y compris à travers fin de mois et année bissextile', () => {
    expect(ajouterJours('2028-02-28', 1)).toBe('2028-02-29')
    expect(ajouterJours('2026-12-31', 1)).toBe('2027-01-01')
    expect(ajouterJours('2026-03-01', -1)).toBe('2026-02-28')
    expect(joursEntre('2026-01-01', '2026-03-01')).toBe(59)
    expect(joursEntre('2026-03-29', '2026-03-28')).toBe(-1)
  })
})

describe('écrit', () => {
  it('compte à rebours seulement avec une date', () => {
    expect(joursAvantEcrit(null, '2026-03-10')).toBeNull()
    expect(joursAvantEcrit('2027-03-04', '2026-03-10')).toBe(359)
  })
})
