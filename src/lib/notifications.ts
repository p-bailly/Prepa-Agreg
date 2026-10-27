// Rappel quotidien par notification locale (facultatif). Programmé par le système Android : il part même appli fermée.
// Dans un navigateur (développement), il n'y a pas de rappel.
import { Capacitor } from '@capacitor/core'
import { LocalNotifications } from '@capacitor/local-notifications'

export const natif = Capacitor.isNativePlatform()
const ID_RAPPEL = 1

/** Programme (ou reprogramme) le rappel quotidien à « HH:MM ». Renvoie faux si la permission est refusée. */
export async function programmerRappel(heure: string): Promise<boolean> {
  if (!natif) return false
  try {
    const perm = await LocalNotifications.requestPermissions()
    if (perm.display !== 'granted') return false
    await LocalNotifications.cancel({ notifications: [{ id: ID_RAPPEL }] })
    const [h, m] = heure.split(':').map(Number)
    await LocalNotifications.schedule({
      notifications: [
        {
          id: ID_RAPPEL,
          title: 'Prépa Agreg',
          body: 'Les cartes du jour vous attendent (10 à 15 minutes).',
          schedule: { on: { hour: h, minute: m }, allowWhileIdle: true },
        },
      ],
    })
    return true
  } catch {
    return false
  }
}

export async function annulerRappel(): Promise<void> {
  if (!natif) return
  try {
    await LocalNotifications.cancel({ notifications: [{ id: ID_RAPPEL }] })
  } catch {
    /* rien à annuler */
  }
}
