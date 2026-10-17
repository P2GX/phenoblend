import { Component, OnInit, signal } from '@angular/core';
import { NavbarComponent } from './navbar/navbar.component';
import { FooterComponent } from 'ng-hpo-uikit';
import { RouterOutlet } from '@angular/router';
import { open as openExternalBrowser } from '@tauri-apps/plugin-shell';
import { getCurrentWindow } from '@tauri-apps/api/window';
import { getVersion } from '@tauri-apps/api/app';


@Component({
  selector: 'app-root',
  standalone: true,
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css', '../styles.scss'],
  imports: [
    FooterComponent,
    NavbarComponent,
    RouterOutlet,

]
})
export class AppComponent implements OnInit {
  protected appVersion = signal<string>("1.0.0"); // fallback

  async ngOnInit() {
      const appWindow = getCurrentWindow();
      try {
        // pulls the version from tauri.conf.json (which syncs with Cargo.toml)
        const version = await getVersion();
        this.appVersion.set(version);
      } catch (error) {
        console.error('Failed to load app version:', error);
      }
      
      // Listen for the event emitted by Rust backend
      await appWindow.listen('close-requested', async () => {
        // 1. Optional: Add an unsaved changes check here if needed
        // 2. Destroy the window cleanly from the frontend
        await appWindow.destroy();
      });
    }
  
  handleHelpNavigation() {
    this.handleExternalNavigation("https://github.com/P2GX/phenoblend");
  }


  private async handleExternalNavigation(url: string): Promise<void> {
    try {
      await openExternalBrowser(url);
    } catch (error) {
      console.warn('Tauri environment missing, falling back to standard web navigation.', error);
      window.open(url, '_blank', 'noopener,noreferrer');
    }
  }
}
