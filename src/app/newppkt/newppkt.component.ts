import { Component, inject, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { HierarchyMapItem, HpoTermMinimal, NotificationService, OntologyMatch, PhenopacketLoaderComponent, PolishedHpoAnnotation } from 'ng-hpo-uikit';
import { from, Observable, of } from 'rxjs';
import { catchError } from 'rxjs/operators';
import { ConfigService } from '../services/config-service';
import { Router } from '@angular/router';
import { AnnotationService } from '../services/annotation-service';
import { HpoDialogWrapperComponent } from '../util/hpotwostep/hpo-dialog-wrapper.component';
import {  HpoTwostepData } from 'ng-hpo-uikit';

/*
 * This component allows the user to enter HPO data by either choosing a phenopacket or by
 * entering a clinical text and using named-entity recognition (NER)/text-mining to
 * get a list of HPO terms. We expect these term to be the clinical manifestations of
 * an individual with multiple genetic diagnoses.
*/
@Component({
  selector: 'app-new-case',
  standalone: true,
  imports: [CommonModule, PhenopacketLoaderComponent, HpoDialogWrapperComponent],
  templateUrl: './newppkt.component.html',
  styleUrls: ['./newppkt.component.scss']
})
export class NewCaseComponent {

  private router = inject(Router);
  readonly isProcessing = signal<boolean>(false);
  readonly statusMessage = signal<string | null>(null);
  readonly errorDetails = signal<string | null>(null);


  private configService = inject(ConfigService);
  private notificationService = inject(NotificationService);
  private annotationService = inject(AnnotationService);
  protected hierarchyCache = signal<Record<string, HierarchyMapItem>>({});

  readonly observedHpoCount: number = this.annotationService.observedHpoCount();

  /**
   * The explicit callback handler passed down to the library loader.
   * Defined as an arrow function to lock contextual execution scope execution bounds to this class.
   */
  readonly processPhenopacketIngest = async (payload: string): Promise<void> => {
    // Reset state states before starting ingest processing
    this.isProcessing.set(true);
    this.statusMessage.set('Parsing phenopacket payload in Rust engine...');
    this.errorDetails.set(null);
    try { 
      await this.configService.ingestPhenopacket(payload);
      const hpoCount = await this.configService.getObservedHpoCount();
        console.log("hpoCount", hpoCount)
      this.notificationService.showSuccess(`Phenopacket input with ${hpoCount} observed HPO teams.`); 
      this.proceedToNextWindow(hpoCount);
    } catch (error) {
      this.notificationService.showError('Ingestion failed.');
       console.log("Ingestion failed", error)
    } finally {
      this.isProcessing.set(false);
    }
  };


  private readonly availableModifiers = (): Promise<HpoTermMinimal[]> => {
    return this.configService.getHpoModifiers();
  };

  performHpoAutocomplete = (query: string): Observable<OntologyMatch[]> => {
    return from(this.configService.performHpoAutocomplete(query)).pipe(
      catchError(err => {
        this.notificationService.showError(String(err));
        return of([]);
      })
    );
  };
   
  fetchHpoHierarchy = (termId: string): Promise<HierarchyMapItem> => {
    const cached = this.hierarchyCache()[termId];
    if (cached) {
      return Promise.resolve(cached);
    }
    return this.configService.getHpoParentAndChildrenTerms(termId).then(data => {
      this.hierarchyCache.update(cache => ({ ...cache, [termId]: data }));
      return data;
    });
  };

  showHpoTwoStepDialog = signal(false);
  protected hpoDialogData?: HpoTwostepData;

  protected openCurationWizard(): void {
    this.hpoDialogData = {
      mineTextProvider: (text: string) => this.configService.mineClinicalText(text),
      autocompleteProvider: (query: string) => this.performHpoAutocomplete(query),
      hierarchyProvider: (termId: string) => this.fetchHpoHierarchy(termId),
      availableModifiers: () => this.availableModifiers()
    };
    this.showHpoTwoStepDialog.set(true);
  }

  onHpoDialogResult(polishedAnnotations?: PolishedHpoAnnotation[]): void {
    this.showHpoTwoStepDialog.set(false);
    if (polishedAnnotations) {
      const observedTerms: PolishedHpoAnnotation[] = polishedAnnotations.filter((annot) => ! annot.excluded);
       const termIds = observedTerms.map(t => t.termId);
         this.configService.addObservedHposFromNER(termIds);
        const n_observed = observedTerms.length;
       if (n_observed > 0) {
        this.proceedToNextWindow(n_observed);
       } else {
          this.notificationService.showError(`Extracted ${polishedAnnotations.length} phenotype annotations but no observed HPOs!`)
       } 
    }
  }


  private proceedToNextWindow(observedHpoCount: number): void {
    this.annotationService.setObservedHpoCount(observedHpoCount)
    this.router.navigate(['/genedisease']);
  }
}