// Charge TOUS les fichiers de content/ (lecture directe, indépendante du chargement de l'appli) et contrôle leur forme,
// leurs liens et chaque formule KaTeX. Voir CLAUDE.md, « Format du contenu ».
import { readdirSync, readFileSync } from 'node:fs'
import { join } from 'node:path'
import katex from 'katex'
import { describe, expect, it } from 'vitest'
import { CONTENU } from './chargement'
import { decouperTexte, formulesDe } from './texte'
import type { Carte, Contenu, Exercice, Theme } from './types'
import { textesDe, validerContenu } from './validation'
import { carte, contenu, exercice } from '../repetition/outils-test'
import type { Exercice as Ex, FichierQcm, IntroLecon } from './types'

/** Exercice de test avec assez d'indices pour sa difficulté. */
const exo = (id: string, c: string, d: number): Ex => ({ ...exercice(id, c, d), indices: d >= 3 ? ['i1', 'i2'] : [] })

const RACINE = join(__dirname, '../../content')

function lireTout(): { contenu: Contenu; erreurs: string[] } {
  const erreurs: string[] = []
  const themes = (JSON.parse(readFileSync(join(RACINE, 'themes.json'), 'utf8')) as { themes: Theme[] }).themes
  const lireDossier = <T,>(dossier: 'cartes' | 'exercices'): T[] => {
    const sortie: T[] = []
    for (const theme of readdirSync(join(RACINE, dossier))) {
      for (const f of readdirSync(join(RACINE, dossier, theme))) {
        // <niveau>.json, ou <niveau>-<chapitre>.json pour un fichier réservé à un chapitre.
        const m = /^(terminale|l1|l2|l3|agreg)(?:-(.+))?\.json$/.exec(f)
        if (!m) {
          erreurs.push(`${dossier}/${theme}/${f} : nom de fichier inattendu (attendu <niveau>.json ou <niveau>-<chapitre>.json)`)
          continue
        }
        const [, niveau, chapitre] = m
        const donnees = JSON.parse(readFileSync(join(RACINE, dossier, theme, f), 'utf8')) as Record<string, T[]>
        for (const x of donnees[dossier] ?? []) {
          sortie.push(x)
          // Le fichier range la carte (ou l'exercice, via sa carte) sous son thème et son niveau.
          if (dossier === 'cartes') {
            const k = x as Carte
            if (k.theme !== theme || k.niveau !== niveau || (chapitre && k.chapitre !== chapitre))
              erreurs.push(`${dossier}/${theme}/${f} : carte ${k.id} mal rangée`)
          }
        }
        if (!Array.isArray(donnees[dossier])) erreurs.push(`${dossier}/${theme}/${f} : clé « ${dossier} » absente`)
      }
    }
    return sortie
  }
  const cartes = lireDossier<Carte>('cartes')
  const exercices = lireDossier<Exercice>('exercices')
  const qcm = readdirSync(join(RACINE, 'qcm'))
    .filter((f) => f.endsWith('.json'))
    .map((f) => {
      const q = JSON.parse(readFileSync(join(RACINE, 'qcm', f), 'utf8')) as FichierQcm
      if (`${q.chapitre}.json` !== f) erreurs.push(`qcm/${f} : le nom du fichier doit être <chapitre>.json`)
      return q
    })
  const lecons = readdirSync(join(RACINE, 'lecons'))
    .filter((f) => f.endsWith('.json'))
    .map((f) => {
      const l = JSON.parse(readFileSync(join(RACINE, 'lecons', f), 'utf8')) as IntroLecon
      if (`${l.chapitre}.json` !== f) erreurs.push(`lecons/${f} : le nom du fichier doit être <chapitre>.json`)
      return l
    })
  const parId = new Map(cartes.map((c) => [c.id, c]))
  for (const e of exercices) {
    const k = parId.get(e.carte)
    if (k && !e.id.startsWith(k.id)) erreurs.push(`exercice ${e.id} : identifiant attendu « ${k.id}-exN »`)
  }
  return { contenu: { themes, cartes, exercices, lecons, qcm }, erreurs }
}

describe('contenu de content/', () => {
  const { contenu: tout, erreurs } = lireTout()

  it('a la forme attendue (identifiants, thèmes, niveaux, prérequis, trous, affirmations, exercices)', () => {
    expect([...erreurs, ...validerContenu(tout)]).toEqual([])
  })

  it('est celui que l’appli embarque', () => {
    expect(CONTENU.cartes.map((c) => c.id).sort()).toEqual(tout.cartes.map((c) => c.id).sort())
    expect(CONTENU.exercices.length).toBe(tout.exercices.length)
    expect(CONTENU.lecons?.length).toBe(tout.lecons?.length)
    expect(CONTENU.qcm?.length).toBe(tout.qcm?.length)
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

  it('n’a pas de formule mal composée que KaTeX accepte', () => {
    // \mathbb8 (valeur substituée par un générateur) ou P^-1 (s'affiche P⁻1) passent KaTeX mais s'affichent mal.
    const fautes = textesDe(tout)
      .flatMap(({ ou, texte }) => formulesDe(texte).map((f) => ({ ou, f: f.t })))
      .filter(({ f }) => /\\mathbb(?!\s*\{?[A-Z])|\^-[0-9A-Za-z]/.test(f))
      .map(({ ou, f }) => `${ou} : ${f}`)
    expect(fautes).toEqual([])
  })

  it('n’a pas de notation LaTeX hors des formules', () => {
    // « X^a − 1 » ou « \mathbb{R} » écrits hors de $…$ s'afficheraient tels quels.
    const fautes = textesDe(tout)
      .flatMap(({ ou, texte }) =>
        decouperTexte(texte).flatMap((p) => (p.k === 'p' ? [p.morceaux] : p.k === 'liste' ? p.items : [])).flat().map((m) => ({ ou, m })),
      )
      .filter(({ m }) => m.k === 'texte' && /[\\^_{}]/.test(m.t))
      .map(({ ou, m }) => `${ou} : ${m.k === 'texte' ? m.t : ''}`)
    expect(fautes).toEqual([])
  })
})

describe('le contrôle du contenu repère les erreurs', () => {
  const base = () => {
    const a = carte('a', 'analyse', 'terminale')
    const b = carte('b', 'analyse', 'l1', ['a'])
    const ex = (c: string) => [exo(`${c}-ex1`, c, 1), exo(`${c}-ex2`, c, 2), exo(`${c}-ex3`, c, 4)]
    return { a, b, ex: [...ex('a'), ...ex('b')] }
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
    const exs = ['c1', 'c2'].flatMap((c) => [exo(`${c}-ex1`, c, 1), exo(`${c}-ex2`, c, 2), exo(`${c}-ex3`, c, 5)])
    expect(validerContenu(contenu([c1, c2], exs)).join()).toMatch(/cycle/)
  })
  it('refuse des questions de QCM mal formées', () => {
    const { a, b, ex } = base()
    const c = contenu([a, b], ex)
    const ch = a.chapitre
    const q = (x: object) => ({ id: `${ch}-q1`, carte: 'a', genre: 'choix', enonce: 'e', propositions: ['1', '2', '3', '4'], bonne: 0, explication: 'x', ...x })
    expect(validerContenu({ ...c, qcm: [{ chapitre: ch, questions: [q({})] }] } as never)).toEqual([])
    const msg = (x: object) => validerContenu({ ...c, qcm: [{ chapitre: ch, questions: [q(x)] }] } as never).join()
    expect(msg({ propositions: ['1', '1', '2', '3'] })).toMatch(/non distinctes/)
    expect(msg({ bonne: 4 })).toMatch(/« bonne » invalide/)
    expect(msg({ carte: 'b' })).toMatch(/absente du chapitre/)
    expect(msg({ explication: '' })).toMatch(/obligatoires/)
    expect(msg({ id: 'autre-q1' })).toMatch(/identifiant attendu/)
    expect(msg({ genre: 'vrai-faux' })).toMatch(/un vrai-faux a « vrai »/)
    expect(validerContenu({ ...c, qcm: [{ chapitre: ch, questions: [{ id: `${ch}-q2`, carte: 'a', genre: 'vrai-faux', enonce: 'e', vrai: false, explication: 'x' }] }] } as never)).toEqual([])
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
    const memes = [exo('a-ex1', 'a', 2), exo('a-ex2', 'a', 2), exo('a-ex3', 'a', 2), ...ex.slice(3)]
    const msg0 = validerContenu(contenu([a, b], memes)).join()
    expect(msg0).toMatch(/difficultés distinctes/)
    expect(msg0).toMatch(/difficulté 4 ou 5/)
    expect(validerContenu(contenu([a, b], ex.slice(1))).join()).toMatch(/3 à 6 exercices/)
    const rep = { ...ex[0], reponse_finale: { question: 'q', propositions: ['1', '2', '2', '3'], bonne: 4 } }
    const msg = validerContenu(contenu([a, b], [rep, ...ex.slice(1)])).join()
    expect(msg).toMatch(/non distinctes/)
    expect(msg).toMatch(/« bonne » invalide/)
    const v = { ...ex[0], reponse_finale: { question: 'q', propositions: ['1', '2', '3', '4'], bonne: 0 } }
    const sansReponse = { ...v, variantes: [{ enonce: 'autre', correction: 'c' }] }
    expect(validerContenu(contenu([a, b], [sansReponse, ...ex.slice(1)])).join()).toMatch(/variante 1 : « reponse_finale »/)
    const memeEnonce = { ...v, variantes: [{ enonce: v.enonce, correction: 'c', reponse_finale: v.reponse_finale }] }
    expect(validerContenu(contenu([a, b], [memeEnonce, ...ex.slice(1)])).join()).toMatch(/même énoncé/)
    const bonne = { ...v, variantes: [{ enonce: 'autre', correction: 'c', reponse_finale: { ...v.reponse_finale, bonne: 2 } }] }
    expect(validerContenu(contenu([a, b], [bonne, ...ex.slice(1)]))).toEqual([])
    const sansIndice = { ...ex[1], difficulte: 3 }
    expect(validerContenu(contenu([a, b], [ex[0], sansIndice, ...ex.slice(2)])).join()).toMatch(/2 indices/)
  })
})
