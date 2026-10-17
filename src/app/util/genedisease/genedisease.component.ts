import { Component, inject, signal, computed, ElementRef, viewChild } from '@angular/core';
import { FormControl, ReactiveFormsModule } from '@angular/forms';
import { toSignal } from '@angular/core/rxjs-interop';
import { startWith, debounceTime, switchMap, from, of, catchError } from 'rxjs';
import { GeneDiseaseAssociation } from '../../models/interfaces';
import { AnnotationService } from '../../services/annotation-service';
import { Router } from '@angular/router';
import { NotificationService } from 'ng-hpo-uikit';
import { ConfigService } from '../../services/config-service';


interface GeneEntry {
  geneSymbol: string;
  associations: GeneDiseaseAssociation[];
  selectedDiseaseIds: ReadonlySet<string>;
}

@Component({
  selector: 'genedisease',
  standalone: true,
  imports: [ReactiveFormsModule],
  templateUrl: './genedisease.component.html',
  styleUrls: ['./genedisease.component.scss']
})
export class GeneDiseaseComponent {
  private readonly annotationService = inject(AnnotationService);
  private configService = inject(ConfigService);
  private readonly router = inject(Router);
  private readonly notificationService = inject(NotificationService);

  protected control = new FormControl<string>('', { nonNullable: true });
  protected isOpen = signal(false);
  protected highlightedIndex = signal(-1);
  protected geneEntries = signal<GeneEntry[]>([]);
  private inputEl = viewChild.required<ElementRef<HTMLInputElement>>('geneInput');

  protected readonly genesWithSelections = computed(() =>
    this.geneEntries().filter(e => e.selectedDiseaseIds.size > 0).length
  );

  

  protected readonly canProceed = computed(() => this.genesWithSelections() >= 2);

  private searchResults = toSignal(
    this.control.valueChanges.pipe(
      startWith(this.control.value),
      debounceTime(300),
      switchMap(query => {
        const trimmed = query.trim();
        if (trimmed.length < 2) {
          return of<GeneDiseaseAssociation[]>([]);
        }
        return from(this.annotationService.autocompleteGeneSymbol(trimmed)).pipe(
          catchError(err => {
            this.notificationService.showError(String(err));
            return of<GeneDiseaseAssociation[]>([]);
          })
        );
      })
    ),
    { initialValue: [] as GeneDiseaseAssociation[] }
  );

  // Distinct gene symbols for the dropdown — backend returns one row per
  // gene+disease pair, so the same gene symbol can appear many times.
  protected geneSymbolOptions = computed(() => {
    const seen = new Set<string>();
    const symbols: string[] = [];
    for (const assoc of this.searchResults()) {
      if (!seen.has(assoc.geneSymbol)) {
        seen.add(assoc.geneSymbol);
        symbols.push(assoc.geneSymbol);
      }
    }
    return symbols;
  });

 

  protected onFocus(): void {
    this.isOpen.set(true);
  }

  protected onInput(): void {
    this.isOpen.set(true);
    this.highlightedIndex.set(-1);
  }

  protected onBlur(): void {
    // Deferred so a (mousedown) selection on an option fires first — see template comment.
    this.isOpen.set(false);
    this.highlightedIndex.set(-1);
  }

  protected onKeyDown(event: KeyboardEvent): void {
    const options = this.geneSymbolOptions();
    if (!options.length) return;

    switch (event.key) {
      case 'ArrowDown':
        event.preventDefault();
        this.isOpen.set(true);
        this.highlightedIndex.update(i => (i + 1) % options.length);
        break;
      case 'ArrowUp':
        event.preventDefault();
        this.isOpen.set(true);
        this.highlightedIndex.update(i => (i - 1 + options.length) % options.length);
        break;
      case 'Enter': {
        const idx = this.highlightedIndex();
        if (idx >= 0 && idx < options.length) {
          event.preventDefault();
          this.selectGene(options[idx]);
        }
        break;
      }
      case 'Escape':
        this.isOpen.set(false);
        this.highlightedIndex.set(-1);
        break;
    }
  }

  protected selectGene(geneSymbol: string): void {
    this.addGene(geneSymbol);
    this.control.setValue('');
    this.isOpen.set(false);
    this.highlightedIndex.set(-1);
    this.inputEl().nativeElement.focus();
  }

  private addGene(geneSymbol: string): void {
    if (this.geneEntries().some(e => e.geneSymbol === geneSymbol)) {
      return; // already added — no duplicates
    }

    const associations = this.searchResults().filter(a => a.geneSymbol === geneSymbol);
    const entry: GeneEntry = {
      geneSymbol,
      associations,
      selectedDiseaseIds: new Set()
    };
    this.geneEntries.update(entries => [...entries, entry]);
  }

  protected toggleDisease(geneSymbol: string, diseaseId: string): void {
    this.geneEntries.update(entries =>
      entries.map(e => {
        if (e.geneSymbol !== geneSymbol) return e;
        const next = new Set(e.selectedDiseaseIds);
        next.has(diseaseId) ? next.delete(diseaseId) : next.add(diseaseId);
        return { ...e, selectedDiseaseIds: next };
      })
    );
  }

  protected isDiseaseSelected(entry: GeneEntry, diseaseId: string): boolean {
    return entry.selectedDiseaseIds.has(diseaseId);
  }

  protected selectAllDiseases(geneSymbol: string): void {
    this.geneEntries.update(entries =>
      entries.map(e => {
        if (e.geneSymbol !== geneSymbol) return e;
        return { ...e, selectedDiseaseIds: new Set(e.associations.map(a => a.diseaseId)) };
      })
    );
  }

  protected clearAllDiseases(geneSymbol: string): void {
    this.geneEntries.update(entries =>
      entries.map(e => (e.geneSymbol === geneSymbol ? { ...e, selectedDiseaseIds: new Set() } : e))
    );
  }

  protected removeGene(geneSymbol: string): void {
    this.geneEntries.update(entries => entries.filter(e => e.geneSymbol !== geneSymbol));
  }

  protected readonly hasAnySelection = computed(() =>
    this.geneEntries().some(e => e.selectedDiseaseIds.size > 0)
  );

  protected goToVisualization(): void {
    const byGene = new Map<string, GeneDiseaseAssociation[]>();

    for (const entry of this.geneEntries()) {
      const chosen : GeneDiseaseAssociation[] = entry.associations.filter(a => entry.selectedDiseaseIds.has(a.diseaseId));
      if (chosen.length > 0) {
        byGene.set(entry.geneSymbol, chosen);
      }
    }

    this.annotationService.setSelectedAssociations(byGene);
    this.configService.initializeNotebook(byGene);
    this.router.navigate(['/visualize']);
  }
}