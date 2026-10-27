import { describe, expect, it } from 'vitest'
import { decouperLigne, decouperTexte, formulesDe, trousDe } from './texte'

describe('découpage du texte', () => {
  it('sépare texte, formules en ligne et gras', () => {
    expect(decouperLigne('Soit $f$ **continue** sur $\\left[a,b\\right]$.')).toEqual([
      { k: 'texte', t: 'Soit ' },
      { k: 'math', t: 'f', bloc: false },
      { k: 'texte', t: ' ' },
      { k: 'texte', t: 'continue', gras: true },
      { k: 'texte', t: ' sur ' },
      { k: 'math', t: '\\left[a,b\\right]', bloc: false },
      { k: 'texte', t: '.' },
    ])
  })
  it('reconnaît les trous et les dollars échappés', () => {
    expect(decouperLigne('f est {{1}} et coûte 3 \\$')).toEqual([
      { k: 'texte', t: 'f est ' },
      { k: 'trou', n: 1 },
      { k: 'texte', t: ' et coûte 3 $' },
    ])
    expect(trousDe('a {{1}} b {{2}} c {{1}}')).toEqual([1, 2, 1])
  })
  it('garde une formule non refermée comme du texte', () => {
    expect(decouperLigne('prix $5')).toEqual([{ k: 'texte', t: 'prix $5' }])
  })
  it('découpe paragraphes, listes et formules centrées', () => {
    const p = decouperTexte('Premier\nsuite.\n\n- un $x$\n- deux\n\n$$\\int_0^1 f$$\n\nFin $$y=1$$ ici.')
    expect(p.map((x) => x.k)).toEqual(['p', 'liste', 'bloc', 'p'])
    expect(formulesDe('Fin $$y=1$$ et $z$.\n\n$$\\sum_k k$$')).toEqual([
      { t: 'y=1', bloc: true },
      { t: 'z', bloc: false },
      { t: '\\sum_k k', bloc: true },
    ])
  })
  it('ne prend pas $ … $ … $$ pour une seule formule centrée', () => {
    expect(decouperTexte('$$a$$ et $$b$$').map((x) => x.k)).toEqual(['p'])
  })
})
