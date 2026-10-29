import { describe, expect, it } from 'vitest'
import { CONTENU, GRAPHE } from './chargement'
import { disposerPyramide, ECART, liensChapitres, validerGraphe, type FichierGraphe } from './graphe'
import { carte, contenu } from '../repetition/outils-test'
import { rangNiveau } from './types'

describe('graphe des chapitres', () => {
  it('content/graphe.json est valide : chapitres connus, liens montants, tout relié, sans cycle', () => {
    expect(validerGraphe(CONTENU, GRAPHE)).toEqual([])
  })

  it('repère les erreurs', () => {
    const c = contenu([carte('a', 'x', 'l1'), carte('b', 'y', 'l2'), carte('c', 'z', 'l2')])
    const g: FichierGraphe = {
      courts: { 'x-l1': 'X', 'y-l2': 'Y' },
      liens: [
        ['y-l2', 'x-l1'],
        ['x-l1', 'x-l1'],
        ['x-l1', 'inconnu'],
      ],
    }
    const err = validerGraphe(c, g).join('\n')
    expect(err).toContain('nom court manquant pour z-l2')
    expect(err).toContain('qui descend de niveau')
    expect(err).toContain('vers lui-même')
    expect(err).toContain('chapitre inconnu')
    expect(err).toContain('z-l2 relié à rien')
    const cycle: FichierGraphe = { courts: { 'x-l1': 'X', 'y-l2': 'Y', 'z-l2': 'Z' }, liens: [['y-l2', 'z-l2'], ['z-l2', 'y-l2'], ['x-l1', 'y-l2']] }
    expect(validerGraphe(c, cycle).join('\n')).toContain('cycle')
  })

  it('ajoute les liens imposés par les prérequis des cartes', () => {
    const c = contenu([carte('a', 'x', 'l1'), carte('b', 'y', 'l2', ['a'])])
    expect(liensChapitres(c, { courts: {}, liens: [] })).toEqual([{ de: 'x-l1', vers: 'y-l2' }])
  })

  it('place chaque chapitre une fois, sans chevauchement, la Terminale en bas et l’agrégation en haut', () => {
    const p = disposerPyramide(CONTENU, liensChapitres(CONTENU, GRAPHE))
    const ids = CONTENU.themes.flatMap((t) => t.chapitres.map((ch) => ch.id))
    expect(p.noeuds.map((n) => n.id).sort()).toEqual([...ids].sort())
    for (const a of p.noeuds)
      for (const b of p.noeuds) if (a !== b) expect(Math.hypot(a.x - b.x, a.y - b.y)).toBeGreaterThanOrEqual(ECART / 2)
    for (const a of p.noeuds)
      for (const b of p.noeuds) if (rangNiveau(a.niveau) < rangNiveau(b.niveau)) expect(a.y).toBeGreaterThan(b.y)
    for (const n of p.noeuds) {
      expect(n.x).toBeGreaterThan(0)
      expect(n.x).toBeLessThan(p.largeur)
      expect(n.y).toBeLessThan(p.hauteur)
    }
  })
})
