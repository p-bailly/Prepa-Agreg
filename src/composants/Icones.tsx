// Icônes dessinées (pas d'emoji) : cadenas du parcours.
export function Cadenas({ taille = 20 }: { taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="2">
      <rect x="5" y="11" width="14" height="10" />
      <path d="M8 11V7a4 4 0 0 1 8 0v4" />
    </svg>
  )
}

/** Pastille d'état : plein (acquis), moitié (en cours), vide (à venir), tiret (masqué). */
export function Pastille({ etat, taille = 22 }: { etat: 'plein' | 'moitie' | 'vide' | 'tiret'; taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" stroke="currentColor" strokeWidth="2">
      {etat === 'tiret' ? (
        <path d="M6 12h12" />
      ) : (
        <>
          <circle cx="12" cy="12" r="9" fill={etat === 'plein' ? 'currentColor' : 'none'} />
          {etat === 'moitie' && <path d="M12 3a9 9 0 0 1 0 18z" fill="currentColor" />}
        </>
      )}
    </svg>
  )
}

/** Case d'un niveau d'exercices : ■ réussi, ◧ fait en partie, □ à faire (dessinées : pas de glyphe à chercher). */
export function CaseNiveau({ etat, taille = 22 }: { etat: 'reussi' | 'partiel' | 'a-faire' | 'aucun'; taille?: number }) {
  return (
    <svg width={taille} height={taille} viewBox="0 0 24 24" aria-hidden="true" stroke="currentColor" strokeWidth="2">
      <rect x="3" y="3" width="18" height="18" fill={etat === 'reussi' ? 'currentColor' : 'none'} />
      {etat === 'partiel' && <rect x="3" y="3" width="9" height="18" fill="currentColor" />}
    </svg>
  )
}
