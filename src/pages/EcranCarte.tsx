// Écran d'une carte : recto selon le mode de révision, puis verso en pile de volets (énoncé, démonstration,
// contre-exemple, exercices). Le volet principal est noté par le pied fixé ; chaque autre volet à revoir, une fois
// ouvert, affiche ses propres boutons de notation.
import { useMemo, useState, type CSSProperties } from 'react'
import { BandeauDoute, Encadre, Volet } from '../composants/Blocs'
import { BoutonsNotation } from '../composants/Notation'
import { CaseNiveau } from '../composants/Icones'
import { Signaler } from '../composants/Signalement'
import { Texte } from '../composants/Texte'
import { cleVolet, NOMS_NIVEAUX, NOMS_VOLETS, voletsDe, type Carte, type Volet as NomVolet } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { exerciceSuggere, niveauxDe, NOMS_ETATS_NIVEAU } from '../exercices/niveaux'
import { melanger } from '../lib/hasard'
import type { Mode } from '../lib/reglages'
import { formaterJours, NOMS_NOTES, type Note } from '../repetition/planification'

interface Props {
  carte: Carte
  /** « n / N » dans la séance ; absent en révision libre. */
  position?: { n: number; total: number }
  /** Volets à revoir, le premier étant le principal. */
  volets: NomVolet[]
  neuve: boolean
  mode: Mode
  surNote: (volet: NomVolet, note: Note, principal: boolean) => void
  surExercice: (id: string) => void
  surQuitter: () => void
  libelleQuitter: string
}

const QUESTIONS: Record<NomVolet, string> = {
  enonce: "Vous avez su l'énoncé ?",
  demo: 'Vous avez su la démonstration ?',
  contre: 'Vous avez su le contre-exemple ?',
}

export function EcranCarte(props: Props) {
  const { carte, position, volets, neuve, mode, surNote, surExercice, surQuitter, libelleQuitter } = props
  const { contenu, progression, horloge, aujourdhui } = useDonnees()
  const [verso, setVerso] = useState(false)
  const principal = volets[0] ?? 'enonce'
  const [ouverts, setOuverts] = useState<Set<NomVolet | 'exercices'>>(() => new Set(['enonce', principal]))
  const [notes, setNotes] = useState<Partial<Record<NomVolet, Note>>>({})
  const theme = contenu.themes.find((t) => t.id === carte.theme)
  const maintenant = horloge()
  const jour = aujourdhui()

  const basculer = (v: NomVolet | 'exercices', o: boolean) =>
    setOuverts((s) => {
      const n = new Set(s)
      if (o) n.add(v)
      else n.delete(v)
      return n
    })

  const noterSecondaire = (v: NomVolet, n: Note) => {
    setNotes((x) => ({ ...x, [v]: n }))
    surNote(v, n, false)
  }

  const style = { '--pied': verso ? '140px' : '100px' } as CSSProperties
  const etatVolet = (v: NomVolet) => progression.volets[cleVolet(carte.id, v)]

  return (
    <div className="appli" style={style}>
      <main className="ecran">
        <div className="barre-haut">
          <button className="retour" onClick={surQuitter}>
            ← {libelleQuitter}
          </button>
          {position ? (
            <span className="avancement">
              <span className="piste" aria-hidden="true">
                <span style={{ width: `${(100 * position.n) / Math.max(1, position.total)}%` }} />
              </span>
              {position.n} / {position.total}
            </span>
          ) : (
            <span className="avancement">Révision libre</span>
          )}
        </div>
        <header className="tete-carte">
          <span className="sur-titre">
            {neuve ? 'Nouvelle carte · ' : ''}
            {theme?.titre ?? carte.theme} · {NOMS_NIVEAUX[carte.niveau]}
          </span>
          <h1 className="titre-carte">
            <Texte texte={carte.titre} enLigne />
          </h1>
        </header>
        <BandeauDoute doute={carte.doute} />

        {!verso && <Recto carte={carte} principal={principal} mode={mode} neuve={neuve} />}

        {verso && (
          <>
            {voletsDe(carte).map((v) => {
              const aRevoir = volets.includes(v)
              const note = notes[v]
              const marque =
                note !== undefined ? `noté : ${NOMS_NOTES[note].toLowerCase()}` : aRevoir ? (etatVolet(v) ? 'à revoir' : 'à voir') : undefined
              return (
                <Volet key={v} titre={NOMS_VOLETS[v]} marque={marque} ouvert={ouverts.has(v)} surBascule={(o) => basculer(v, o)}>
                  <ContenuVolet carte={carte} volet={v} />
                  {aRevoir && v !== principal && (
                    <div className="notation-volet">
                      {note === undefined ? (
                        <>
                          <p className="question-pied">{QUESTIONS[v]}</p>
                          <BoutonsNotation etat={etatVolet(v)} maintenant={maintenant} jour={jour} surNote={(n) => noterSecondaire(v, n)} />
                        </>
                      ) : (
                        <p className="petit doux">
                          Prochaine révision de ce volet :{' '}
                          {note === 0 ? 'dans la séance' : `dans ${formaterJours(etatVolet(v)?.intervalle ?? 0)}`}.
                        </p>
                      )}
                    </div>
                  )}
                </Volet>
              )
            })}
            <Volet titre="Exercices" ouvert={ouverts.has('exercices')} surBascule={(o) => basculer('exercices', o)}>
              <VoletExercices carte={carte} surExercice={surExercice} />
            </Volet>
          </>
        )}
        <Signaler cible={carte.id} genre="carte" ou={verso ? undefined : 'recto'} />
      </main>

      <div className="pied">
        <div className="pied-interieur">
          {!verso ? (
            <button className="bouton-principal" onClick={() => setVerso(true)}>
              Révéler
            </button>
          ) : (
            <>
              <p className="question-pied">{QUESTIONS[principal]}</p>
              <BoutonsNotation
                etat={etatVolet(principal)}
                maintenant={maintenant}
                jour={jour}
                surNote={(n) => surNote(principal, n, true)}
              />
            </>
          )}
        </div>
      </div>
    </div>
  )
}

/** Contenu d'un volet noté : l'énoncé dans son encadré (et l'exemple), la démonstration, le contre-exemple. */
function ContenuVolet({ carte, volet }: { carte: Carte; volet: NomVolet }) {
  if (volet === 'enonce')
    return (
      <>
        <Encadre genre={carte.type}>
          <Texte texte={carte.enonce} />
        </Encadre>
        {carte.exemple && (
          <Encadre genre="exemple">
            <Texte texte={carte.exemple} />
          </Encadre>
        )}
      </>
    )
  if (volet === 'demo')
    return (
      <div style={{ paddingTop: 12 }}>
        <Texte texte={carte.demo ?? ''} />
      </div>
    )
  return (
    <Encadre genre="remarque" nom="Hypothèses nécessaires">
      <Texte texte={carte.contre_exemple ?? ''} />
    </Encadre>
  )
}

/** Recto : la consigne, et l'exercice propre au mode (hypothèses à cocher, trous, vrai/faux). */
function Recto({ carte, principal, mode, neuve }: { carte: Carte; principal: NomVolet; mode: Mode; neuve: boolean }) {
  if (principal === 'demo' || principal === 'contre')
    return (
      <>
        <p className="consigne">
          {principal === 'demo'
            ? 'Démontrer ce résultat sur papier ou de tête, puis révéler la démonstration.'
            : 'Pour chaque hypothèse importante, trouver un contre-exemple qui montre qu’elle est nécessaire, puis révéler.'}
        </p>
        <Encadre genre={carte.type}>
          <Texte texte={carte.enonce} />
        </Encadre>
      </>
    )
  if (mode === 'hypotheses') return <ModeHypotheses carte={carte} />
  if (mode === 'trous') return <ModeTrous carte={carte} />
  if (mode === 'vraifaux') return <ModeVraiFaux carte={carte} />
  return (
    <>
      <p className="consigne">
        <Texte texte={carte.recto} enLigne />
      </p>
      {neuve && (
        <p className="petit doux">
          Nouvelle carte : si vous ne la connaissez pas encore, révélez-la et lisez-la attentivement, hypothèses comprises.
        </p>
      )}
    </>
  )
}

/** Verdict d'une ligne : icône et texte, jamais la couleur seule. */
function Verdict({ juste, texte }: { juste: boolean; texte: string }) {
  return (
    <span className={`verdict ${juste ? 'juste' : 'fausse'}`}>
      {juste ? '✓ Juste' : '✗ Faux'} : {texte}
    </span>
  )
}

function ModeHypotheses({ carte }: { carte: Carte }) {
  const lignes = useMemo(
    () =>
      melanger([
        ...(carte.hypotheses ?? []).map((t) => ({ t, vraie: true })),
        ...(carte.hypotheses_pieges ?? []).map((t) => ({ t, vraie: false })),
      ]),
    [carte],
  )
  const [coches, setCoches] = useState<Set<number>>(new Set())
  const [verifie, setVerifie] = useState(false)
  const justes = lignes.filter((l, i) => coches.has(i) === l.vraie).length
  return (
    <>
      <p className="consigne">Cocher les hypothèses du résultat « <Texte texte={carte.titre} enLigne /> », et elles seules.</p>
      <ul className="liste-choix">
        {lignes.map((l, i) => {
          const coche = coches.has(i)
          const juste = coche === l.vraie
          return (
            <li key={i}>
              <button
                className={`case${verifie ? (juste ? ' juste' : ' fausse') : ''}`}
                role="checkbox"
                aria-checked={coche}
                disabled={verifie}
                onClick={() =>
                  setCoches((s) => {
                    const n = new Set(s)
                    if (n.has(i)) n.delete(i)
                    else n.add(i)
                    return n
                  })
                }
              >
                <span className="coche" aria-hidden="true">
                  {coche ? '✓' : ''}
                </span>
                <span className="contenu">
                  <Texte texte={l.t} enLigne />
                  {verifie && (
                    <Verdict
                      juste={juste}
                      texte={l.vraie ? 'c’est une hypothèse du résultat' : 'ce n’est pas une hypothèse du résultat'}
                    />
                  )}
                </span>
              </button>
            </li>
          )
        })}
      </ul>
      {!verifie ? (
        <button className="bouton" onClick={() => setVerifie(true)}>
          Vérifier
        </button>
      ) : (
        <p>
          <strong>
            {justes} / {lignes.length}
          </strong>{' '}
          lignes justes. Révélez l’énoncé pour le relire en entier.
        </p>
      )}
    </>
  )
}

function ModeTrous({ carte }: { carte: Carte }) {
  const trous = carte.trous ?? []
  const options = useMemo(
    () => trous.map((t) => melanger([{ t: t.reponse, bonne: true }, ...t.leurres.map((l) => ({ t: l, bonne: false }))])),
    [trous],
  )
  const [choix, setChoix] = useState<(number | null)[]>(() => trous.map(() => null))
  const [verifie, setVerifie] = useState(false)
  const choisi = (n: number) => {
    const c = choix[n - 1]
    return c === null || c === undefined ? null : options[n - 1][c]
  }
  return (
    <>
      <p className="consigne">Compléter l’énoncé : choisir la bonne proposition pour chaque trou.</p>
      <Encadre genre={carte.type}>
        <Texte
          texte={carte.enonce_trous ?? ''}
          trou={(n) => {
            const c = choisi(n)
            const etat = verifie && c ? (c.bonne ? ' juste' : ' fausse') : ''
            return (
              <span className={`trou${c ? ' rempli' : ''}${etat}`}>
                {c ? <Texte texte={c.t} enLigne /> : n}
                {verifie && c && <span aria-label={c.bonne ? 'juste' : 'faux'}>{c.bonne ? ' ✓' : ' ✗'}</span>}
              </span>
            )
          }}
        />
      </Encadre>
      {options.map((ops, k) => (
        <div key={k} style={{ marginBottom: 14 }}>
          <span className="etiquette">Trou {k + 1}</span>
          <ul className="liste-choix">
            {ops.map((o, i) => {
              const pris = choix[k] === i
              const classe = verifie ? (o.bonne ? ' juste' : pris ? ' fausse' : '') : pris ? ' choisie' : ''
              return (
                <li key={i}>
                  <button
                    className={`case${classe}`}
                    aria-pressed={pris}
                    disabled={verifie}
                    onClick={() => setChoix((c) => c.map((x, j) => (j === k ? i : x)))}
                  >
                    <span className="coche" aria-hidden="true">
                      {pris ? '●' : ''}
                    </span>
                    <span className="contenu">
                      <Texte texte={o.t} enLigne />
                      {verifie && o.bonne && <Verdict juste texte="la bonne réponse" />}
                      {verifie && pris && !o.bonne && <Verdict juste={false} texte="ce leurre rend l’énoncé faux" />}
                    </span>
                  </button>
                </li>
              )
            })}
          </ul>
        </div>
      ))}
      {!verifie && (
        <button className="bouton" disabled={choix.some((c) => c === null)} onClick={() => setVerifie(true)}>
          Vérifier
        </button>
      )}
    </>
  )
}

function ModeVraiFaux({ carte }: { carte: Carte }) {
  const [a] = useState(() => melanger(carte.affirmations ?? [])[0])
  const [reponse, setReponse] = useState<boolean | null>(null)
  const [voir, setVoir] = useState(false)
  if (!a) return null
  return (
    <>
      <p className="consigne">Vrai ou faux ? Si c’est faux, trouver un contre-exemple.</p>
      <Encadre genre="remarque" nom="Affirmation">
        <Texte texte={a.texte} />
      </Encadre>
      {reponse === null ? (
        <div className="deux-colonnes">
          <button className="bouton" onClick={() => setReponse(true)}>
            Vrai
          </button>
          <button className="bouton" onClick={() => setReponse(false)}>
            Faux
          </button>
        </div>
      ) : (
        <>
          <p>
            <Verdict juste={reponse === a.vrai} texte={a.vrai ? 'l’affirmation est vraie.' : 'l’affirmation est fausse.'} />
          </p>
          {!a.vrai && !voir && (
            <>
              <p className="petit doux">Retrouvez un contre-exemple sur papier, puis comparez.</p>
              <button className="bouton" onClick={() => setVoir(true)}>
                Voir le contre-exemple
              </button>
            </>
          )}
          {!a.vrai && voir && a.contre_exemple && (
            <Encadre genre="exemple" nom="Contre-exemple">
              <Texte texte={a.contre_exemple} />
            </Encadre>
          )}
          {a.vrai && a.justification && (
            <Encadre genre="remarque" nom="Pourquoi">
              <Texte texte={a.justification} />
            </Encadre>
          )}
        </>
      )}
    </>
  )
}

/** Volet Exercices : cinq cases de difficulté, l'exercice suggéré, la liste. */
export function VoletExercices({ carte, surExercice }: { carte: Carte; surExercice: (id: string) => void }) {
  const { exercicesParCarte, progression } = useDonnees()
  const exercices = exercicesParCarte.get(carte.id) ?? []
  const niveaux = niveauxDe(exercices, progression.exercices)
  const suggere = exerciceSuggere(exercices, progression.exercices)
  if (!exercices.length) return <p className="doux">Pas encore d’exercice pour cette notion.</p>
  return (
    <div style={{ paddingTop: 12 }}>
      <div className="niveaux-exercices">
        {niveaux.map((n) => (
          <button
            key={n.difficulte}
            disabled={n.etat === 'aucun'}
            className={suggere && suggere.difficulte === n.difficulte ? 'propose' : undefined}
            aria-label={`Difficulté ${n.difficulte} : ${NOMS_ETATS_NIVEAU[n.etat]}`}
            onClick={() => {
              const e = n.exercices.find((x) => x.id === suggere?.id) ?? n.exercices[0]
              if (e) surExercice(e.id)
            }}
          >
            <CaseNiveau etat={n.etat} />
            <span>{n.difficulte}</span>
          </button>
        ))}
      </div>
      <p className="petit doux">■ réussi · ◧ fait en partie · □ à faire · grisé : pas d’exercice à ce niveau</p>
      {suggere ? (
        <button className="bouton" onClick={() => surExercice(suggere.id)}>
          Exercice proposé (difficulté {suggere.difficulte}) : <Texte texte={suggere.titre} enLigne />
        </button>
      ) : (
        <p>Tous les exercices de cette notion sont réussis.</p>
      )}
      <ul className="liste-cartes" style={{ marginTop: 12 }}>
        {exercices.map((e) => (
          <li key={e.id}>
            <button className="bouton compact" style={{ width: '100%', justifyContent: 'flex-start', textAlign: 'left' }} onClick={() => surExercice(e.id)}>
              <span className="doux">{e.difficulte}/5</span>
              <span>
                <Texte texte={e.titre} enLigne />
              </span>
            </button>
          </li>
        ))}
      </ul>
    </div>
  )
}
