// QCM d'un chapitre : dix questions tirées au hasard (vrai ou faux, énoncé à compléter, question d'exercice), une à la
// fois ; après chaque réponse, la bonne est montrée avec sa justification. À la fin, le score est enregistré.
import { useMemo, useState, type CSSProperties } from 'react'
import { chapitresNumerotes } from '../chapitres/plan'
import { CONSIGNES, questionsChapitre, tirerQcm } from '../chapitres/qcm'
import { Volet } from '../composants/Blocs'
import { Texte } from '../composants/Texte'
import { useDonnees } from '../etat/Donnees'
import { enregistrerQcm } from '../lib/progression'

export function EcranQcm({ chapitre, surFin }: { chapitre: string; surFin: () => void }) {
  const [tirage, setTirage] = useState(0)
  return <Deroule key={tirage} chapitre={chapitre} surFin={surFin} recommencer={() => setTirage(tirage + 1)} />
}

function Deroule({ chapitre, surFin, recommencer }: { chapitre: string; surFin: () => void; recommencer: () => void }) {
  const { contenu, cartesParId, changerProgression, horloge, aujourdhui } = useDonnees()
  const entree = useMemo(() => chapitresNumerotes(contenu).find((x) => x.chapitre.id === chapitre), [contenu, chapitre])
  const questions = useMemo(() => tirerQcm(questionsChapitre(contenu, chapitre)), [contenu, chapitre])
  const [rang, setRang] = useState(0)
  const [choix, setChoix] = useState<(number | null)[]>(() => questions.map(() => null))
  const fini = rang >= questions.length
  const bonnes = questions.filter((q, i) => choix[i] === q.bonne).length
  const q = questions[rang]
  const pris = q ? choix[rang] : null

  const suivante = () => {
    if (rang === questions.length - 1)
      changerProgression((p) => enregistrerQcm(p, chapitre, { date: horloge(), jour: aujourdhui(), bonnes, total: questions.length }))
    setRang(rang + 1)
    window.scrollTo(0, 0)
  }

  const style = { '--pied': '100px' } as CSSProperties
  return (
    <div className="appli" style={style}>
      <main className="ecran">
        <div className="barre-haut">
          <button className="retour" onClick={surFin}>
            ← {fini ? 'Chapitre' : 'Quitter'}
          </button>
          <span className="avancement">
            <span className="piste" aria-hidden="true">
              <span style={{ width: `${(100 * Math.min(rang, questions.length)) / Math.max(1, questions.length)}%` }} />
            </span>
            {fini ? 'Résultat' : `${rang + 1} / ${questions.length}`}
          </span>
        </div>
        <span className="sur-titre-ecran">
          QCM · Chapitre {entree?.numero} · {entree?.chapitre.titre}
        </span>

        {q && (
          <>
            <span className="etiquette">
              {CONSIGNES[q.genre]} · notion « <Texte texte={cartesParId.get(q.carte)?.titre ?? ''} enLigne /> »
            </span>
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
                  ? 'Très bien : chapitre maîtrisé à ce QCM.'
                  : 'Revoyez les notions concernées, puis retentez : les questions changent à chaque tirage.'}
            </p>
            {questions.some((x, i) => choix[i] !== x.bonne) && (
              <>
                <span className="etiquette">À revoir</span>
                <ul className="liste-revoir">
                  {[...new Set(questions.filter((x, i) => choix[i] !== x.bonne).map((x) => x.carte))].map((id) => (
                    <li key={id}>
                      <Texte texte={cartesParId.get(id)?.titre ?? id} enLigne />
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
                Retour au chapitre
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
