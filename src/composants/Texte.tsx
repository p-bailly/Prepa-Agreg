// Affiche un texte du contenu : paragraphes, listes, gras, formules KaTeX (en ligne et centrées, ces dernières défilant
// dans leur cadre si elles sont trop larges), trous {{n}} des énoncés à trous. Voir CLAUDE.md, « Texte et formules ».
import katex from 'katex'
import { Fragment, memo, type ReactNode } from 'react'
import { decouperLigne, decouperTexte, type Morceau } from '../contenu/texte'

const cache = new Map<string, string>()

function formule(t: string, bloc: boolean): string {
  const cle = `${bloc ? 'B' : 'L'}${t}`
  let html = cache.get(cle)
  if (html === undefined) {
    // throwOnError : false à l'affichage (le test du contenu, lui, exige que tout passe).
    html = katex.renderToString(t, { throwOnError: false, displayMode: bloc, strict: 'ignore', output: 'html' })
    cache.set(cle, html)
  }
  return html
}

export function Formule({ t, bloc = false }: { t: string; bloc?: boolean }) {
  if (bloc) return <span className="formule-bloc" dangerouslySetInnerHTML={{ __html: formule(t, true) }} />
  return <span dangerouslySetInnerHTML={{ __html: formule(t, false) }} />
}

function Morceaux({ morceaux, trou }: { morceaux: Morceau[]; trou?: (n: number) => ReactNode }) {
  return (
    <>
      {morceaux.map((m, i) => {
        if (m.k === 'math') return <Formule key={i} t={m.t} bloc={m.bloc} />
        if (m.k === 'trou') return <Fragment key={i}>{trou ? trou(m.n) : `[${m.n}]`}</Fragment>
        return m.gras ? <strong key={i}>{m.t}</strong> : <Fragment key={i}>{m.t}</Fragment>
      })}
    </>
  )
}

interface Props {
  texte: string
  /** Rendu d'un trou {{n}}. */
  trou?: (n: number) => ReactNode
  /** Une seule ligne, sans paragraphe (boutons, listes de choix). */
  enLigne?: boolean
}

export const Texte = memo(function Texte({ texte, trou, enLigne }: Props) {
  if (enLigne) return <Morceaux morceaux={decouperLigne(texte.replace(/\s*\n\s*/g, ' '))} trou={trou} />
  return (
    <div className="texte">
      {decouperTexte(texte).map((p, i) => {
        if (p.k === 'bloc') return <Formule key={i} t={p.t} bloc />
        if (p.k === 'liste')
          return (
            <ul key={i}>
              {p.items.map((it, j) => (
                <li key={j}>
                  <Morceaux morceaux={it} trou={trou} />
                </li>
              ))}
            </ul>
          )
        return (
          <p key={i}>
            <Morceaux morceaux={p.morceaux} trou={trou} />
          </p>
        )
      })}
    </div>
  )
})
