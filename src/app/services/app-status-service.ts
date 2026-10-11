import { Injectable, signal, inject, NgZone, computed } from '@angular/core';
import { listen } from '@tauri-apps/api/event';
import { invoke } from '@tauri-apps/api/core';
import { NotificationService, OntologyLoadEvent } from 'ng-hpo-uikit';
import { WritableSignal } from '@angular/core';
import { InitializationStatusDto } from '../models/status_dto'; 
import { initialNotebook, Notebook } from '../models/notebook';


@Injectable({ providedIn: 'root' })
export class AppStatusService {

  private ngZone = inject(NgZone);
  private notificationService = inject(NotificationService);

  hpoLoaded = signal<boolean>(false);
  hpoaLoaded = signal<boolean>(false);
  g2dLoaded = signal<boolean>(false);
  hpoVersion = signal<string>('');
  hpoaVersion = signal<string>('');
  g2dVersion = signal<string>('');
  hpoJsonPath = signal<string>('');
  nHpoTerms = signal<number>(0);
  nHpoaDisease = signal<number>(0);
  hpoLoading = signal<boolean>(false);
  hpoaLoading = signal<boolean>(false);
  g2dLoading = signal<boolean>(false);
  nGenes = signal<number>(0);
  biocuratorOrcid = signal<string>('');
  errorMessage = signal<string>('');
  hasError = computed(() => !!this.errorMessage());

  readonly notebook = signal<Notebook | null>(null);

  constructor() {
    this.setupListeners().then(() => {
      this.fetchPreloadedStatus();
      this.registerNotebook();
    });
  }


 private async setupListeners() {
    // HPO
    await this.registerOntologyListener({
      channel: 'hpo-load-event',
      loadingSignal: this.hpoLoading,
      loadedSignal: this.hpoLoaded,
      versionSignal: this.hpoVersion,
      countSignal: this.nHpoTerms,
      errorContext: 'HPO'
    });

    // HPOA 
    await this.registerOntologyListener({
      channel: 'hpoa-load-event',
      loadingSignal: this.hpoaLoading,
      loadedSignal: this.hpoaLoaded,
      versionSignal: this.hpoaVersion,
      countSignal: this.nHpoaDisease,
      errorContext: 'HPOA'
    });

    // Gene-Disease associations
    await this.registerOntologyListener({
      channel: 'g2d-load-event',
      loadingSignal: this.g2dLoading,
      loadedSignal: this.g2dLoaded,
      versionSignal: this.g2dVersion,
      countSignal: this.nGenes,
      errorContext: 'gene-to-disease'
    });
  }




private async registerNotebook(): Promise<void> {
   await listen<Notebook>('notebook-updated', (event) => {
      this.ngZone.run(() => this.notebook.set(event.payload));
    });
    try {
      const initial = initialNotebook();
      if (this.notebook() === null) {
        this.ngZone.run(() => this.notebook.set(initial));
      }
    } catch (err) {
      console.error('Failed to fetch initial notebook:', err);
    }
}


/**
 * Registers a generic Tauri event listener for ontology loading (HPO or MAxO)
 */
private async registerOntologyListener(config: {
  channel: string;
  loadingSignal: WritableSignal<boolean>;
  loadedSignal: WritableSignal<boolean>;
  versionSignal: WritableSignal<string>;
  countSignal: WritableSignal<number>;
  errorContext: string;
}) {
 
  await listen(config.channel, (event) => {
    const { status, payload } = event.payload as OntologyLoadEvent;

    this.ngZone.run(() => {
      switch (status) {
        case 'loading':
          config.loadingSignal.set(true);
          this.errorMessage.set('');
          break;

        case 'success':
          config.loadingSignal.set(false);
          config.loadedSignal.set(true);
          const versionInfo = payload?.version || 'n/a';
          const totalTerms = payload?.termCount ?? 0;
          config.countSignal.set(totalTerms);
          config.versionSignal.set(versionInfo);
          break;

        case 'error':
          config.loadingSignal.set(false);
          const errorMsg = typeof payload === 'object' ? payload?.errorMessage : payload;
          this.errorMessage.set(`[${config.errorContext}] ${errorMsg || 'Failed to parse'}`);
          this.notificationService.showError(this.errorMessage());
          break;

        case 'cancel':
          config.loadingSignal.set(false);
          break;
      }
    });
  });
}

  private async fetchPreloadedStatus() {
    try {
      const status = await invoke<InitializationStatusDto>('check_initialization_status');
      
      this.ngZone.run(() => {
        this.hpoLoaded.set(status.hpo_loaded);
        this.nHpoTerms.set(status.hpo_terms);

        this.hpoaLoaded.set(status.hpoa_loaded);
        this.nHpoaDisease.set(status.hpoa_diseases);

        this.g2dLoaded.set(status.g2d_loaded);
        this.nGenes.set(status.n_genes);
      });
    } catch (err: any) {
      console.error('Failed to resolve initial backend pre-load state context:', err);
    }
  }
    

}