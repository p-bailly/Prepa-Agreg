// Écriture du fichier d'export : dans l'APK, fichier temporaire puis feuille de partage Android (Drive, mail, Fichiers…) ;
// dans un navigateur (développement), téléchargement.
import { Capacitor } from '@capacitor/core'

export async function partagerFichier(nom: string, contenu: string): Promise<void> {
  if (Capacitor.isNativePlatform()) {
    const { Filesystem, Directory, Encoding } = await import('@capacitor/filesystem')
    const { Share } = await import('@capacitor/share')
    const ecrit = await Filesystem.writeFile({ path: nom, data: contenu, directory: Directory.Cache, encoding: Encoding.UTF8 })
    await Share.share({ title: nom, url: ecrit.uri, dialogTitle: 'Enregistrer la sauvegarde' })
    return
  }
  const url = URL.createObjectURL(new Blob([contenu], { type: 'application/json' }))
  const a = document.createElement('a')
  a.href = url
  a.download = nom
  a.click()
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}
