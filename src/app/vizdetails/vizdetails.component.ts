import { ChangeDetectionStrategy, Component, computed, input, output } from '@angular/core';
import { PercentPipe } from '@angular/common';
import { UpsetPlotData } from 'ngx-phenoprofile';
import { GeneMatch, Notebook } from '../models/notebook';


export type VisualizationType = 'overlap' | 'upset' | 'spread';
export type DetailsExportFormat = 'json' | 'html';

/** Short definitions shown with the details. Move to a constants file when it grows. */
const MATCH_TYPE_HELP: { type: GeneMatch['matchType']; text: string }[] = [
  { type: 'exact', text: 'The gene’s disease is annotated to exactly this term.' },
  {
    type: 'descendant',
    text: 'The disease is annotated to a more specific term (a descendant) of this term.',
  },
  { type: 'partial', text: 'No direct match; the score reflects the similarity of a shared ancestor.' },
  { type: 'noMatch', text: 'Only the root phenotype term is shared.' },
];

/**
 * The explanation and the underlying data for the plot that is currently shown.
 *
 * Presentational only: it renders what the notebook contains and asks the parent to
 * produce downloads (so file saving can use the Tauri dialog in one place).
 */
@Component({
  selector: 'app-viz-details',
  standalone: true,
  imports: [PercentPipe],
  changeDetection: ChangeDetectionStrategy.OnPush,
  templateUrl: './vizdetails.component.html',
  styleUrls: [ './vizdetails.component.scss' ],
})
export class VizDetailsComponent {
  view = input.required<VisualizationType>();
  notebook = input.required<Notebook>();
  /** The parent saves the file (Tauri save dialog) so this component stays presentational. */
  exportRequested = output<DetailsExportFormat>();

  readonly matchTypeHelp = MATCH_TYPE_HELP;

  readonly overlap = computed(() => this.notebook().analysis?.overlap);
  readonly upset = computed(() => this.notebook().analysis?.upset as UpsetPlotData | undefined);
  readonly spread = computed(() => this.notebook().analysis?.spread);

  /** "NF1+PTPN11" for the genes at `indices` of `upset.rows`. */
  genesOf(indices: number[], upset: UpsetPlotData): string {
    return indices.map((i) => upset.rows[i].geneSymbol).join('+');
  }

  segmentLabel(indices: number[], upset: UpsetPlotData): string {
    return indices.length === 0 ? 'Not explained' : this.genesOf(indices, upset);
  }

  /** One overlap cell as text; a gene absent from the row has no match recorded. */
  describeMatch(m: GeneMatch | undefined): string {
    if (!m) {
      return 'no match';
    }
    const base = `${m.matchType} (${m.score})`;
    return m.matchType === 'noMatch' ? base : `${base}: ${m.matchedTerm.termId} ${m.matchedTerm.termLabel}`;
  }
}
