import {
  Component,
  ViewChild,
  signal,
  computed,
  inject,
  OnInit,
} from '@angular/core';
import * as d3 from 'd3';
import { ConfigService } from '../services/config-service';
import { AnnotationService } from '../services/annotation-service';
import { DuoEnrichmentSummary, EnrichmentResult, GeneDiseaseAssociation } from '../models/interfaces';
import { NotificationService } from 'ng-hpo-uikit';
import {
  UpsetPlotComponent,
  SpreadPlotComponent,
  OverlapPlotComponent,
  SpreadPlotPayload,
  OverlapPlotPayload,
  UpsetPlotPayload,
} from 'ngx-phenoprofile';

//  3 view modes
type VisualizationType = 'overlap' | 'upset' | 'spread';
type ExportFormat = 'svg' | 'pdf';

@Component({
  selector: 'app-presence-visualizer',
  standalone: true,
  imports: [OverlapPlotComponent, UpsetPlotComponent, SpreadPlotComponent],
  templateUrl: './visualize.component.html',
  styleUrls: ['./visualize.component.scss'],
})
export class PhenotypeProfileVisualizerComponent implements OnInit {
  @ViewChild('overlapComponent') private childMatrix!: OverlapPlotComponent;
  @ViewChild('upsetComponent') private childUpset!: UpsetPlotComponent;
  @ViewChild('spreadComponent') private childSpread!: SpreadPlotComponent;

  private configService = inject(ConfigService);
  protected annotationService = inject(AnnotationService);
  private notificationService = inject(NotificationService);

  readonly matrixData = signal<OverlapPlotPayload | null>(null);
  readonly upsetData = signal<UpsetPlotPayload | null>(null);
  readonly spreadData = signal<SpreadPlotPayload | null>(null);
  readonly isLoading = signal<boolean>(false);
  readonly activeView = signal<VisualizationType>('overlap');

  enrichmentResults = signal<DuoEnrichmentSummary[] | null>(null);
  isEnrichmentLoading = signal(false);
  enrichmentError = signal<string | null>(null);

  hasEnrichmentResults = computed(() => (this.enrichmentResults()?.length ?? 0) > 0);


  // Computed title depending on which view is selected
  readonly currentTitle = computed(() => {
    switch (this.activeView()) {
      case 'overlap':
        return 'Phenotypic Overlap Plot';
      case 'upset':
        return 'Phenotypic Upset Plot';
      case 'spread':
        return 'Phenotypic Spread Plot';
    }
  });

  // Handler for drop-down change event - change view type
  onViewChange(event: Event): void {
    const target = event.target as HTMLSelectElement;
    this.activeView.set(target.value as VisualizationType);
  }

  ngOnInit(): void {
    this.reloadMatrixData();
  }

  async reloadMatrixData(): Promise<void> {
    this.isLoading.set(true);
    try {
      const selectedAnnotations: GeneDiseaseAssociation[] =
        this.annotationService.allSelectedAssociations();
      if (selectedAnnotations.length < 2) {
        const errMsg = `At least two gene/disease pairs required to perform analysis but only ${selectedAnnotations.length} available.`;
        this.notificationService.showError(errMsg);
        return;
      }
      const annotationMap = selectedAnnotations.reduce((acc, item) => {
        const key = item.geneSymbol;
        if (!acc.has(key)) {
          acc.set(key, []);
        }
        acc.get(key)!.push(item);
        return acc;
      }, new Map<string, GeneDiseaseAssociation[]>());
 
      
      await this.configService.initializeNotebook(annotationMap);
      const overlapResult =
        await this.configService.getOverlapPlotData();
      const upsetResult =
        await this.configService.getUpsetPlotPayload();
      const spreadResult =
        await this.configService.getSpreadPlotPayload();
      this.upsetData.set(upsetResult);
      this.matrixData.set(overlapResult);
      this.spreadData.set(spreadResult);
    } catch (err) {
      this.notificationService.showError(
        `Failed fetching matrix values: ${err}.`,
      );
    } finally {
      this.isLoading.set(false);
    }
  }

  async exportMatrix(format: ExportFormat): Promise<void> {
    const view = this.activeView();
    if (view !== 'overlap' && view !== 'upset' && view !== 'spread') return;

    const targetElement = this.getChartContainerElement(view);
    if (!targetElement) {
      this.notificationService.showError('Visualization element not found.');
      return;
    }

    const svgElement = d3
      .select(targetElement)
      .select('svg')
      .node() as SVGElement | null;
    if (!svgElement) {
      this.notificationService.showError(
        'No SVG content found for the active view.',
      );
      return;
    }

    const svgClone = svgElement.cloneNode(true) as SVGElement;
    svgClone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svgClone.setAttribute('version', '1.1');

    const serializer = new XMLSerializer();
    const svgString =
      '<?xml version="1.0" standalone="no"?>\n' +
      serializer.serializeToString(svgClone);

    const dateStamp = new Date().toISOString().split('T')[0];
    const plotType = this.activeView(); // overlap, spread, etc.
    const baseName = this.annotationService.defaultFileName(plotType);
    const filenameBase = `${baseName}_${dateStamp}`;

    try {
      switch (format) {
        case 'svg': {
          const saved = await this.configService.saveSvgFile(
            svgString,
            `${filenameBase}.svg`,
          );
          if (saved)
            this.notificationService.showSuccess('SVG exported successfully.');
          break;
        }
        case 'pdf': {
          const saved = await this.configService.exportSvgToPdf(
            svgString,
            `${filenameBase}.pdf`,
          );
          if (saved)
            this.notificationService.showSuccess('PDF exported successfully.');
          break;
        }
      }
    } catch (error) {
      console.error(`${format.toUpperCase()} export failed:`, error);
      this.notificationService.showError(
        `Failed to export as ${format.toUpperCase()}.`,
      );
    }
  }

  private getChartContainerElement(
    view: VisualizationType,
  ): HTMLElement | null {
    switch (view) {
      case 'overlap':
        return this.childMatrix?.chartContainerRef?.nativeElement ?? null;
      case 'upset':
        return this.childUpset?.chartContainerRef?.nativeElement ?? null;
      case 'spread':
        return this.childSpread?.chartContainerRef?.nativeElement ?? null;
    }
  }

  async runDuoEnrichment(): Promise<void> {
    this.isEnrichmentLoading.set(true);
    this.enrichmentError.set(null);
    try {
        const recordData = Object.fromEntries(
        this.annotationService.selectedAssociations(),
      );
      const results = await this.configService.analyzeDuoEnrichment(recordData, 10_000);
      this.enrichmentResults.set(results);
    } catch (err) {
      console.error('Duo enrichment analysis failed:', err);
      this.enrichmentError.set(
        err instanceof Error ? err.message : 'Enrichment analysis failed'
      );
      this.enrichmentResults.set(null);
    } finally {
      this.isEnrichmentLoading.set(false);
    }
  }

  /** Height percentages for a simple inline pmf bar chart, 0-100 scale. */
  pmfBarHeights(result: EnrichmentResult): number[] {
    const max = Math.max(...result.pmf, 1e-9); // avoid divide-by-zero on an all-zero pmf
    return result.pmf.map((p) => (p / max) * 100);
  }

  formatPValue(p: number): string {
    return p < 0.001 ? p.toExponential(2) : p.toFixed(3);
  }

  isDegenerate(result: DuoEnrichmentSummary): boolean {
    return result.unweighted?.degenerate ?? result.frequency?.degenerate ?? false;
  }
}
