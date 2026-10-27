// Charge TOUS les fichiers de content/ (lecture directe, indépendante du chargement de l'appli) et contrôle leur forme,
// leurs liens et chaque formule KaTeX. Voir CLAUDE.md, « Format du contenu ».
import { readdirSync, readFileSync } from 'node:fs'
import { join } from 'node:path'
import katex from 'katex'
import { describe, expect, it } from 'vitest'
import { CONTENU } from './chargement'
import { formulesDe } from './texte'
import { NIVEAUX, type Carte, type Contenu, type Exercice, type Theme } from './types'
import { textesDe, validerContenu } from './validation'
import { carte, contenu, exercice } from '../repetition/outils-test'

const RACINE = join(__dirname, '../../content')

function lireTout(): { contenu: Contenu; erreurs: string[] } {
  const erreurs: string[] = []
  const themes = (JSON.parse(readFileSync(join(RACINE, 'themes.json'), 'utf8')) as { themes: Theme[] }).themes
  const lireDossier = <T,>(dossier: 'cartes' | 'exercices'): T[] => {
    const sortie: T[] = []
    for (const theme of readdirSync(join(RACINE, dossier))) {
      for (const f of readdirSync(join(RACINE, dossier, theme))) {
        const niveau = f.replace(/\.json$/, '')
        if (!f.endsWith('.json') || !NIVEAUX.includes(niveau as never)) {
          erreurs.push(`${dossier}/${theme}/${f} : nom de fichier inattendu (attendu <niveau>.json)`)
          continue
        }
        const donnees = JSON.parse(readFileSync(join(RACINE, dossier, theme, f), 'utf8')) as Record<string, T[]>
        for (const x of donnees[dossier] ?? []) {
          sortie.push(x)
          // Le fichier range la carte (ou l'exercice, via sa carte) sous son thème et son niveau.
          if (dossier === 'cartes') {
            const k = x as Carte
            if (k.theme !== theme || k.niveau !== niveau) erreurs.push(`${dossier}/${theme}/${f} : carte ${k.id} mal rangée`)
          }
        }
        if (!Array.isArray(donnees[dossier])) erreurs.push(`${dossier}/${theme}/${f} : clé « ${dossier} » absente`)
      }
    }
    return sortie
  }
  const cartes = lireDossier<Carte>('cartes')
  const exercices = lireDossier<Exercice>('exercices')
  const parId = new Map(cartes.map((c) => [c.id, c]))
  for (const e of exercices) {
    const k = parId.get(e.carte)
    if (k && !e.id.startsWith(k.id)) erreurs.push(`exercice ${e.id} : identifiant attendu « ${k.id}-exN »`)
  }
  return { contenu: { themes, cartes, exercices }, erreurs }
}

describe('contenu de content/', () => {
  const { contenu: tout, erreurs } = lireTout()

  it('a la forme attendue (identifiants, thèmes, niveaux, prérequis, trous, affirmations, exercices)', () => {
    expect([...erreurs, ...validerContenu(tout)]).toEqual([])
  })

  it('est celui que l’appli embarque', () => {
    expect(CONTENU.cartes.map((c) => c.id).sort()).toEqual(tout.cartes.map((c) => c.id).sort())
    expect(CONTENU.exercices.length).toBe(tout.exercices.length)
  })

  it('a chaque formule rendue par KaTeX sans erreur', () => {
    const fautes: string[] = []
    for (const { ou, texte } of textesDe(tout)) {
      for (const f of formulesDe(texte)) {
        try {
          katex.renderToString(f.t, { throwOnError: true, displayMode: f.bloc, strict: 'ignore' })
        } catch (e) {
          fautes.push(`${ou} : ${(e as Error).message}`)
        }
      }
      // Un $ isolé (formule non refermée) s'afficherait tel quel.
      const sansEchappes = texte.replace(/\\\$/g, '')
      if ((sansEchappes.match(/\$/g) ?? []).length % 2) fautes.push(`${ou} : nombre impair de $`)
    }
    expect(fautes).toEqual([])
  })

  it('n’écrit pas les intervalles à l’anglaise', () => {
    const fautes = textesDe(tout)
      .flatMap(({ ou, texte }) => formulesDe(texte).map((f) => ({ ou, f: f.t })))
      .filter(({ f }) => /\(\s*-?\s*\\infty|\\infty\s*\)/.test(f))
      .map(({ ou, f }) => `${ou} : ${f}`)
    expect(fautes).toEqual([])
  })
})

describe('le contrôle du contenu repère les erreurs', () => {
  const base = () => {
    const a = carte('a', 'analyse', 'terminale')
    const b = carte('b', 'analyse', 'l1', ['a'])
    return { a, b, ex: [exercice('a-ex1', 'a', 1), exercice('a-ex2', 'a', 2), exercice('b-ex1', 'b', 1), exercice('b-ex2', 'b', 2)] }
  }
  it('accepte un contenu correct', () => {
    const { a, b, ex } = base()
    expect(validerContenu(contenu([a, b], ex))).toEqual([])
  })
  it('refuse prérequis inconnu, de niveau supérieur, en cycle', () => {
    const { a, b, ex } = base()
    expect(validerContenu(contenu([{ ...a, prerequis: ['zz'] }, b], ex)).join()).toMatch(/prérequis inconnu/)
    expect(validerContenu(contenu([{ ...a, prerequis: ['b'] }, b], ex)).join()).toMatch(/niveau supérieur/)
    const c1 = carte('c1', 'analyse', 'terminale', ['c2'])
    const c2 = carte('c2', 'analyse', 'terminale', ['c1'])
    const exs = [exercice('c1-ex1', 'c1', 1), exercice('c1-ex2', 'c1', 2), exercice('c2-ex1', 'c2', 1), exercice('c2-ex2', 'c2', 2)]
    expect(validerContenu(contenu([c1, c2], exs)).join()).toMatch(/cycle/)
  })
  it('refuse trous mal formés et leurres en double', () => {
    const { a, b, ex } = base()
    const mal = { ...a, enonce_trous: 'f est {{2}}', trous: [{ reponse: 'continue', leurres: ['x', 'y', 'z'] }] }
    expect(validerContenu(contenu([mal, b], ex)).join()).toMatch(/mal numérotés/)
    const double = { ...a, enonce_trous: 'f est {{1}}', trous: [{ reponse: 'x', leurres: ['x', 'y', 'z'] }] }
    expect(validerContenu(contenu([double, b], ex)).join()).toMatch(/non distincts/)
    const dansFormule = { ...a, enonce_trous: '$f = {{1}}$', trous: [{ reponse: 'w', leurres: ['x', 'y', 'z'] }] }
    expect(validerContenu(contenu([dansFormule, b], ex)).join()).toMatch(/à l'intérieur d'une formule/)
  })
  it('refuse vrai/faux sans affirmation fausse, exercices trop peu variés, réponse mal formée', () => {
    const { a, b, ex } = base()
    const vf = { ...a, affirmations: [{ texte: 'x', vrai: true }] }
    expect(validerContenu(contenu([vf, b], ex)).join()).toMatch(/au moins une affirmation fausse/)
    const memes = [exercice('a-ex1', 'a', 2), exercice('a-ex2', 'a', 2), ex[2], ex[3]]
    expect(validerContenu(contenu([a, b], memes)).join()).toMatch(/difficultés distinctes/)
    const rep = { ...ex[0], reponse_finale: { question: 'q', propositions: ['1', '2', '2', '3'], bonne: 4 } }
    const msg = validerContenu(contenu([a, b], [rep, ...ex.slice(1)])).join()
    expect(msg).toMatch(/non distinctes/)
    expect(msg).toMatch(/« bonne » invalide/)
    const sansIndice = { ...ex[1], difficulte: 3 }
    expect(validerContenu(contenu([a, b], [ex[0], sansIndice, ex[2], ex[3]])).join()).toMatch(/2 indices/)
  })
})
