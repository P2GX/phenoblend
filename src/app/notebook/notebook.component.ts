// notebook.component.ts
import { Component, computed, inject, Signal} from '@angular/core';
import { DatePipe, PercentPipe,JsonPipe } from '@angular/common';
import { AppStatusService } from '../services/app-status-service';
import { GeneDiseaseEntity, GeneMatch, OverlapData, PhenoBlendAnalysis, SimpleOntologyTerm, TermMatches, UpsetData } from '../models/notebook';

@Component({
  selector: 'app-notebook',
  standalone: true,
  imports: [DatePipe, JsonPipe, PercentPipe],
  templateUrl: './notebook.component.html',
  styleUrls: ['./notebook.component.scss'],
})
export class NotebookComponent {
  private appService = inject(AppStatusService);
  readonly notebook = this.appService.notebook;
  // Assumes Unix seconds; drop the * 1000 if your timestamp is already milliseconds
  readonly updatedAt = computed(() => {
    const nb = this.notebook();
    return nb ? new Date(nb.timestamp * 1000) : null;
  });
  readonly analysis: Signal<PhenoBlendAnalysis | null> = computed(() => {
    return this.notebook()?.analysis ?? null;
  });
  readonly overlap: Signal<OverlapData|null> = computed(() => {
    return this.analysis()?.overlap ?? null;
  })
  readonly spread = computed(() => this.analysis()?.spread ?? null);

  readonly geneToDiseaseEntityList: Signal<GeneDiseaseEntity[]> = computed<GeneDiseaseEntity[]>(
      () => this.notebook()?.geneDiseaseList ?? []
  );


  readonly geneToHpoListMap = computed<Record<string, SimpleOntologyTerm[]>>(
    () => this.analysis()?.geneAnnotations ?? {}
  );

   // Look up the GeneMatch for a specific entity/gene symbol under this HPO term column
  protected getEntityMatchForTerm(termCol: TermMatches, entity: string): GeneMatch | undefined {
    return termCol.matches[entity];
  }

  // Returns CSS class based on the MatchType
  protected getMatchCellClass(match: GeneMatch | undefined): string {
    if (!match) return 'cell-unmatched';
    
    switch (match.matchType) {
      case 'exact':
        return 'cell-exact';
      case 'descendant':
        return 'cell-descendant';
      case 'partial':
        return 'cell-partial';
      case 'noMatch':
      default:
        return 'cell-unmatched';
    }
  }

  // Generates rich content/text to display inside the cell
  protected getCellDisplayLabel(match: GeneMatch | undefined): string {
    if (!match || match.matchType === 'noMatch') return 'X';
    
    switch (match.matchType) {
      case 'exact':
        return 'Exact';
      case 'descendant':
        // Show type and the specific matched term label/id
        return `Descendant: ${match.matchedTerm.termLabel || match.matchedTerm.termId}`;
      case 'partial':
        // Show type and the partial term label/id
        return `Partial: ${match.matchedTerm.termLabel || match.matchedTerm.termId} (${match.score.toFixed(2)})`;
      default:
        return 'Match';
    }
  }

  protected getMatchTooltip(match: GeneMatch | undefined, entity: string, termLabel: string): string {
    if (!match || match.matchType === 'noMatch') return `${entity} has no match for ${termLabel}`;
    return `${entity} matched via ${match.matchType} (Score: ${match.score})`;
  }

  // Extracts unique HPO terms that are covered/overlapped by at least one gene
  protected getAccountedTerms(upsetData: UpsetData): SimpleOntologyTerm[] {
    const termMap = new Map<string, SimpleOntologyTerm>();
    
    for (const gene of upsetData.genes) {
      for (const match of gene.probandOverlapped) {
        // We track the original proband term that got overlapped
        termMap.set(match.term.termId, match.term);
      }
    }
    
    return Array.from(termMap.values());
  }
 

}