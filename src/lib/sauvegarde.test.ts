import { describe, expect, it } from 'vitest'
import { ajouterSignalement, noterVolet, PROGRESSION_VIDE } from './progression'
import { REGLAGES_DEFAUT } from './reglages'
import { fabriquerSauvegarde, lireSauvegarde } from './sauvegarde'

const T = Date.parse('2026-03-01T09:00:00Z')

describe('export et import de la progression', () => {
  it('fait l’aller-retour, signalements compris', () => {
    let p = noterVolet(PROGRESSION_VIDE, 'c', 'enonce', 2, T, '2026-03-01', true)
    p = ajouterSignalement(p, { date: T, cible: 'c', genre: 'carte', commentaire: 'hypothèse en trop ?' })
    const r = { ...REGLAGES_DEFAUT, plafondNeuves: 7, dateEcrit: '2027-03-04' }
    const texte = JSON.stringify(fabriquerSauvegarde(p, r, T, '1.0.0'))
    const lu = lireSauvegarde(texte)
    expect(lu.progression).toEqual(p)
    expect(lu.reglages).toEqual(r)
  })
  it('refuse un fichier qui n’est pas une sauvegarde', () => {
    expect(() => lireSauvegarde('pas du json')).toThrow(/JSON/)
    expect(() => lireSauvegarde('{"a":1}')).toThrow(/sauvegarde/)
  })
  it('répare des réglages invalides au lieu de les recopier', () => {
    const texte = JSON.stringify({
      format: 'prepa-agreg/progression-v1',
      progression: {},
      reglages: { plafondNeuves: -5, niveauDepart: 'master', modesActifs: ['inconnu'], dateEcrit: 'demain' },
    })
    const { reglages } = lireSauvegarde(texte)
    expect(reglages.plafondNeuves).toBe(0)
    expect(reglages.niveauDepart).toBe('terminale')
    expect(reglages.modesActifs).toEqual(['tete'])
    expect(reglages.dateEcrit).toBeNull()
  })
})
