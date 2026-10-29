// QCM d'un chapitre (dix questions de tous genres) ou d'un niveau entier (vingt questions à choix, ou vingt vrai-faux,
// réparties entre ses chapitres), une question à la fois ; après chaque réponse, la bonne est montrée avec sa
// justification. À la fin, le score est enregistré sous la clé de la cible.
import { useMemo, useState, type CSSProperties } from 'react'
import { chapitresNumerotes } from '../chapitres/plan'
import { CONSIGNES, cleQcm, tirerCible, type CibleQcm } from '../chapitres/qcm'
import { BandeauDoute, Volet } from '../composants/Blocs'
import { Texte } from '../composants/Texte'
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { enregistrerQcm } from '../lib/progression'

export function EcranQcm({ cible, surFin }: { cible: CibleQcm; surFin: () => void }) {
  const [tirage, setTirage] = useState(0)
  return <Deroule key={tirage} cible={cible} surFin={surFin} recommencer={() => setTirage(tirage + 1)} />
}

function Deroule({ cible, surFin, recommencer }: { cible: CibleQcm; surFin: () => void; recommencer: () => void }) {
  const { contenu, cartesParId, changerProgression, horloge, aujourdhui } = useDonnees()
  const numeros = useMemo(() => new Map(chapitresNumerotes(contenu).map((x) => [x.chapitre.id, x])), [contenu])
  const questions = useMemo(() => tirerCible(contenu, cible), [contenu, cible])
  const parChapitre = 'chapitre' in cible
  const entree = parChapitre ? numeros.get(cible.chapitre) : undefined
  const titre = parChapitre
    ? `QCM · Chapitre ${entree?.numero} · ${entree?.chapitre.titre}`
    : `${cible.sorte === 'qcm' ? 'QCM' : 'Vrai ou faux'} · ${NOMS_NIVEAUX[cible.niveau]}, tout le niveau`
  /** Notion d'une question, précédée du numéro de son chapitre pour un QCM de niveau. */
  const notion = (chapitre: string, carte: string) =>
    `${parChapitre ? '' : `chapitre ${numeros.get(chapitre)?.numero} · `}notion « ${cartesParId.get(carte)?.titre ?? carte} »`
  const [rang, setRang] = useState(0)
  const [choix, setChoix] = useState<(number | null)[]>(() => questions.map(() => null))
  const fini = rang >= questions.length
  const bonnes = questions.filter((q, i) => choix[i] === q.bonne).length
  const q = questions[rang]
  const pris = q ? choix[rang] : null

  const suivante = () => {
    if (rang === questions.length - 1)
      changerProgression((p) => enregistrerQcm(p, cleQcm(cible), { date: horloge(), jour: aujourdhui(), bonnes, total: questions.length }))
    setRang(rang + 1)
    window.scrollTo(0, 0)
  }

  const style = { '--pied': '100px' } as CSSProperties
  return (
    <div className="appli" style={style}>
      <main className="ecran">
        <div className="barre-haut">
          <button className="retour" onClick={surFin}>
            ← {fini ? 'Retour' : 'Quitter'}
          </button>
          <span className="avancement">
            <span className="piste" aria-hidden="true">
              <span style={{ width: `${(100 * Math.min(rang, questions.length)) / Math.max(1, questions.length)}%` }} />
            </span>
            {fini ? 'Résultat' : `${rang + 1} / ${questions.length}`}
          </span>
        </div>
        <span className="sur-titre-ecran">{titre}</span>

        {q && (
          <>
            <span className="etiquette">
              <Texte texte={`${CONSIGNES[q.genre]} · ${notion(q.chapitre, q.carte)}`} enLigne />
            </span>
            <BandeauDoute doute={q.doute} />
            {q.contexte && (
              <div className="bloc">
                <Texte texte={q.contexte} />
              </div>
            )}
            <div className="question-qcm">
              <Texte texte={q.enonce} />
            </div>
            <ul className="liste-choix">
              {q.propositions.map((t, i) => {
                const classe = pris === null ? '' : i === q.bonne ? ' juste' : pris === i ? ' fausse' : ''
                return (
                  <li key={i}>
                    <button
                      className={`case${classe}`}
                      aria-pressed={pris === i}
                      disabled={pris !== null}
                      onClick={() => setChoix(choix.map((c, j) => (j === rang ? i : c)))}
                    >
                      <span className="coche" aria-hidden="true">
                        {pris === i ? '●' : ''}
                      </span>
                      <span className="contenu">
                        <Texte texte={t} enLigne />
                        {pris !== null && i === q.bonne && <span className="verdict juste">✓ Bonne réponse</span>}
                        {pris === i && i !== q.bonne && <span className="verdict fausse">✗ Mauvaise réponse</span>}
                      </span>
                    </button>
                  </li>
                )
              })}
            </ul>
            {pris !== null && q.explication && (
              <Volet titre={q.genre === 'calcul' ? 'Correction' : 'Explication'} ouvert={q.genre !== 'calcul' ? true : undefined}>
                <div style={{ paddingTop: 12 }}>
                  <Texte texte={q.explication} />
                </div>
              </Volet>
            )}
          </>
        )}

        {fini && (
          <div className="bloc">
            <p className="grand-nombre">
              {bonnes} / {questions.length}
            </p>
            <p className="doux">
              {bonnes === questions.length
                ? 'Sans faute.'
                : bonnes / questions.length >= 0.8
                  ? `Très bien : ${parChapitre ? 'chapitre maîtrisé' : 'niveau bien en main'} à ce QCM.`
                  : 'Revoyez les notions concernées, puis retentez : les questions changent à chaque tirage.'}
            </p>
            {questions.some((x, i) => choix[i] !== x.bonne) && (
              <>
                <span className="etiquette">À revoir</span>
                <ul className="liste-revoir">
                  {[...new Map(questions.filter((x, i) => choix[i] !== x.bonne).map((x) => [x.carte, x])).values()].map((x) => (
                    <li key={x.carte}>
                      <Texte texte={parChapitre ? (cartesParId.get(x.carte)?.titre ?? x.carte) : `Chapitre ${numeros.get(x.chapitre)?.numero} : ${cartesParId.get(x.carte)?.titre ?? x.carte}`} enLigne />
                    </li>
                  ))}
                </ul>
              </>
            )}
          </div>
        )}
      </main>
      <div className="pied">
        <div className="pied-interieur">
          {fini ? (
            <div className="rangee-boutons">
              <button className="bouton" onClick={recommencer}>
                Nouveau tirage
              </button>
              <button className="bouton" onClick={surFin}>
                {parChapitre ? 'Retour au chapitre' : 'Retour'}
              </button>
            </div>
          ) : (
            <button className="bouton-principal" disabled={pris === null} onClick={suivante}>
              {rang === questions.length - 1 ? 'Voir le résultat' : 'Question suivante'}
            </button>
          )}
        </div>
      </div>
    </div>
  )
}
