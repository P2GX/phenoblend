import { Component, ElementRef, Input, OnChanges, SimpleChanges, ViewChild, signal } from '@angular/core';
import { DecimalPipe} from '@angular/common'
import * as d3 from 'd3';
import { OverlapPlotPayload, Severity } from '../models/phenoprofile_dto';






@Component({
  selector: 'app-overlap-plot',
  standalone: true,
  template: `
    <div class="matrix-container" style="position: relative;">
      <div #chartContainer></div>
      @if (tooltipData(); as data) {
        <div class="matrix-tooltip" [style.left.px]="data.x" [style.top.px]="data.y">
          <strong>{{ data.gene }}</strong> ↔ <strong>{{ data.hpoId }}</strong><br>
          <span class="label">{{ data.hpoName }}</span><br>
          <span class="score">Match Value: {{ data.score | number:'1.2-2' }}</span>
          @if (data.severity) {
            <br><span class="sev">Severity: {{ data.severity }}</span>
          }
        </div>
      }
    </div>
  `,
  styles: [`
    .matrix-container { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }
    .matrix-tooltip {
      position: absolute;
      padding: 8px 12px;
      background: rgba(33, 37, 41, 0.95);
      color: #fff;
      border-radius: 4px;
      font-size: 12px;
      pointer-events: none;
      box-shadow: 0 2px 8px rgba(0,0,0,0.15);
      z-index: 10;
      transform: translate(-50%, -110%);
    }
    .matrix-tooltip .label { color: #adb5bd; font-style: italic; }
    .matrix-tooltip .score { color: #66C2A5; font-weight: bold; }
    .matrix-tooltip .sev { color: #c4b5fd; font-weight: bold; }
  `],
  imports: [DecimalPipe]
})
export class OverlapPlotComponent implements OnChanges {
  @ViewChild('chartContainer', { static: true }) chartContainerRef!: ElementRef<HTMLElement>;

  @Input({ required: true }) data!: OverlapPlotPayload;

 // readonly tooltipData = signal<{ x: number; y: number; gene: string; hpoId: string; hpoName: string; score: number } | null>(null);
  readonly tooltipData = signal<{
      x: number; 
      y: number; 
      gene: string; 
      hpoId: string; 
      hpoName: string; 
      score: number; 
      severity: Severity | null;
  } | null>(null);

  // Match your exact Matplotlib configuration profile colors
  private readonly MATCH_COLOR = '#66C2A5';
  private readonly SIM_COLOR = '#ffd92f';
  private readonly MISMATCH_COLOR = '#ffffff';
  

private readonly SEVERITY_STYLE: Record<Severity, { fill: string; text: string; abbr: string }> = {
  borderline: { fill: '#ede9fe', text: '#4c1d95', abbr: 'Bord.' },
  mild:       { fill: '#c4b5fd', text: '#4c1d95', abbr: 'Mild' },
  moderate:   { fill: '#8b5cf6', text: '#ffffff', abbr: 'Mod.' },
  severe:     { fill: '#5b21b6', text: '#ffffff', abbr: 'Sev.' },
  profound:   { fill: '#2e1065', text: '#ffffff', abbr: 'Prof.' },
};
private readonly HEAVY = new Set<Severity>(['severe', 'profound']);


  ngOnChanges(changes: SimpleChanges): void {
    if (changes['data'] && this.data) {
      this.renderMatrixChart();
    }
  }

  private static readonly LABEL_FONT = '500 11px sans-serif'; // must match the text style below
  private static readonly LABEL_FONT_PX = 11;

  /** The text that is drawn for a term; used both to draw and to measure. */
  private labelText(name: string): string {
    return name.length > 42 ? name.substring(0, 40) + '...' : name;
  }

  /** Room needed by labels rotated -45° that start at the top edge of a column. */
  private rotatedLabelExtent(labels: string[]): { vertical: number; horizontal: number } {
    const ctx = document.createElement('canvas').getContext('2d')!;
    ctx.font = OverlapPlotComponent.LABEL_FONT;
    const longest = Math.max(0, ...labels.map((t) => ctx.measureText(t).width));
    const along = Math.ceil(longest * Math.SQRT1_2); // extent up and to the right
    return {
      vertical: along + Math.ceil(OverlapPlotComponent.LABEL_FONT_PX * Math.SQRT1_2) + 10,
      horizontal: along,
    };
  }

  private renderMatrixChart(): void {
    // Clear out any stale layouts if rendering updates sequentially
    const container = this.chartContainerRef.nativeElement;
    container.innerHTML = '';

    const genes = this.data.entities;
    const cols = this.data.columns;

    // Configurable scaling step bounds (similar to your 0.5 inches space calculations)
    const cellSize = 30;
    const extent = this.rotatedLabelExtent(cols.map((c) => this.labelText(c.hpoName)));
    const margin = {
      top: Math.max(120, extent.vertical),
      right: Math.max(180, extent.horizontal - cellSize / 2 + 10), // the last label slants to the right
      bottom: 70,
      left: 180,
    };
    const matrixWidth = cols.length * cellSize;
    const width = matrixWidth + margin.left + margin.right;
    const height = genes.length * cellSize + margin.top + margin.bottom;

    // 1. Initialize the root SVG canvas viewport
    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height);

    const mainGroup = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // 2. Set up discrete positional linear scales
    const xScale = d3.scaleBand()
      .domain(cols.map(r => r.hpoId))
      .range([0, cols.length * cellSize])
      .padding(0.05);

    const yScale = d3.scaleBand()
      .domain(genes)
      .range([0, genes.length * cellSize])
      .padding(0.05);

    // 3. Custom color interpolator 
    const colorScale = (score: number): string => {
      if (score >= 1.0) return this.MATCH_COLOR;
      if (score <= 0.0) return this.MISMATCH_COLOR;
      
      // Discrete steps for partial patches (Adjust thresholds as needed!)
      if (score < 0.33) return '#fff2ae'; // Very pale yellow
      if (score < 0.66) return '#fede58'; // Mid yellow
      return this.SIM_COLOR;              // Full yellow (#ffd92f)
    };

    // 4. Render X Axis Headers (HPO Terms turned -45 deg)
    const xLabels = mainGroup.append('g')
      .selectAll('.x-label')
      .data(cols)
      .enter()
      .append('g')
      .attr('transform', d => `translate(${xScale(d.hpoId)! + xScale.bandwidth() / 2}, 0)`);

    xLabels.append('text')
      .attr('transform', 'rotate(-45)')
      .attr('dx', '-0.6em')
      .attr('dy', '-0.5em')
      .style('text-anchor', 'start')
      .style('font-size', '11px')
      .style('font-weight', '500')
      .text(d => d.hpoName.length > 42 ? d.hpoName.substring(0, 40) + '...' : d.hpoName);

    // 5. Render Y Axis Headers (Gene Symbols)
    mainGroup.append('g')
      .selectAll('.y-label')
      .data(genes)
      .enter()
      .append('text')
      .attr('x', -12)
      .attr('y', d => yScale(d)! + yScale.bandwidth() / 2)
      .attr('dy', '0.35em')
      .style('text-anchor', 'end')
      .style('font-size', '12px')
      .style('font-weight', 'bold')
      .text(d => d);

    // 6. Draw the Presence Circle Grid Nodes
    cols.forEach((row) => {
      genes.forEach((gene) => {
        const score = row.scores[gene] ?? 0.0;
        const cx = xScale(row.hpoId)! + xScale.bandwidth() / 2;
        const cy = yScale(gene)! + yScale.bandwidth() / 2;

        const cellGroup = mainGroup.append('g')
          .attr('class', 'cell-node');

        const heavy = !!row.severity && score >= 1.0 && this.HEAVY.has(row.severity);
        const baseStroke = heavy ? 2.5 : 1;
        const strokeColor = heavy ? this.SEVERITY_STYLE[row.severity!].fill : '#333333';

        // Render an invisible tracking rectangle wrapper to make hovering easy for small dots
        cellGroup.append('rect')
          .attr('x', xScale(row.hpoId)!)
          .attr('y', yScale(gene)!)
          .attr('width', xScale.bandwidth())
          .attr('height', yScale.bandwidth())
          .style('fill', 'transparent')
          .style('cursor', 'pointer')
          .on('mouseover', (event: MouseEvent) => {
            const [x, y] = d3.pointer(event, container);
            // Highlight row/column intersections dynamically by shifting opacity
            d3.select((event.currentTarget as Element).parentNode as Element).select('circle')
              .style('stroke-width', `${baseStroke+1}px`)
              .attr('r', (xScale.bandwidth() / 2) - 1);

            this.tooltipData.set({
              x, y, gene,
              hpoId: row.hpoId, hpoName: row.hpoName,  score,
              severity: row.severity ?? null,
            });
          })
          .on('mousemove', (event: MouseEvent) => {
             const [x, y] = d3.pointer(event, container);
            this.tooltipData.update(current => current ? { ...current, x, y } : null);
          })
          .on('mouseleave', (event: { currentTarget: { parentNode: any; }; }) => {
            d3.select(event.currentTarget.parentNode).select('circle')
              .style('stroke-width', `${baseStroke}px`)
              .attr('r', (xScale.bandwidth() / 2) - 3);
            this.tooltipData.set(null);
          });

        // The actual visual circle indicator node
        cellGroup.append('circle')
          .attr('cx', cx)
          .attr('cy', cy)
          .attr('r', (xScale.bandwidth() / 2) - 3)
          .style('fill', colorScale(score))
          .style('stroke', strokeColor)
          .style('stroke-width', `${baseStroke}px`)
          .style('pointer-events', 'none'); // Let pointer events fall directly through to the tracking rect
      });
    });

  if (cols.some(c => c.severity)) {
  const trackH = 16;
  const track = mainGroup.append('g')
    .attr('transform', `translate(0, ${genes.length * cellSize + 10})`);

  track.append('text')
    .attr('x', -12).attr('y', trackH / 2).attr('dy', '0.35em')
    .style('text-anchor', 'end').style('font-size', '12px').style('font-weight', 'bold')
    .text('Severity');

  cols.filter(c => c.severity).forEach(c => {
    const st = this.SEVERITY_STYLE[c.severity!];
    const g = track.append('g').attr('transform', `translate(${xScale(c.hpoId)!}, 0)`);
    g.append('rect')
      .attr('width', xScale.bandwidth()).attr('height', trackH).attr('rx', 3)
      .style('fill', st.fill).style('stroke', '#333333').style('stroke-width', '1px');
    g.append('text')
      .attr('x', xScale.bandwidth() / 2).attr('y', trackH / 2).attr('dy', '0.35em')
      .style('text-anchor', 'middle').style('font-size', '9px').style('font-weight', '600')
      .style('fill', st.text).text(st.abbr);
    g.append('title').text(`${c.hpoName}: ${c.severity}`);
  });
}
    // ==========================================
// 7. Render Compact Colorbar Legend on the Right
// ==========================================

const legendData = [
  { label: 'match',   color: this.MATCH_COLOR,   showLabel: true },
  { label: '',        color: '#ffd92f',          showLabel: false }, 
  { label: 'partial', color: '#fede58',          showLabel: true }, 
  { label: '',        color: '#fff2ae',          showLabel: false }, 
  { label: 'none',    color: this.MISMATCH_COLOR, showLabel: true }
];

const legendX = matrixWidth + 40; 
const totalLegendHeight = cellSize * 2; // Hard cap the height to exactly 2 rows (~60px)
const legendItemHeight = totalLegendHeight / legendData.length; // ~12px per block
const legendWidth = 15; // Slightly slimmer to match the compact height

// Create a container group for the legend
const legendGroup = mainGroup.append('g')
  .attr('transform', `translate(${legendX}, 0)`);

// Draw the color rectangles
legendGroup.selectAll('.legend-rect')
  .data(legendData)
  .enter()
  .append('rect')
  .attr('x', 0)
  .attr('y', (_, i) => i * legendItemHeight)
  .attr('width', legendWidth)
  .attr('height', legendItemHeight)
  .style('fill', d => d.color)
  .style('stroke', '#333333')
  .style('stroke-width', '1px');

// Add tick marks on the right side of every block boundary
legendGroup.selectAll('.legend-tick')
  .data(d3.range(1, legendData.length)) // 1 to 4 ensures top and bottom outer borders aren't doubled
  .enter()
  .append('line')
  .attr('x1', 0)
  .attr('x2', legendWidth)
  .attr('y1', d => d * legendItemHeight)
  .attr('y2', d => d * legendItemHeight)
  .style('stroke', '#333333')
  .style('stroke-width', '1px');

// Add text labels only for match, partial, and miss
legendGroup.selectAll('.legend-label')
  .data(legendData)
  .enter()
  .filter(d => d.showLabel) // Only render text for the main milestones
  .append('text')
  .attr('x', legendWidth + 10)
  .attr('y', (_, i, nodes) => {
    // Elegant centering: 'match' aligns with top block, 'miss' with bottom block, 'partial' dead center
    const label = d3.select(nodes[i]).datum() as any;
    if (label.label === 'match') return (legendItemHeight / 2);
    if (label.label === 'none') return totalLegendHeight - (legendItemHeight / 2);
    return totalLegendHeight / 2; // 'partial' centered right in the middle of the bar
  })
  .attr('dy', '0.35em')
  .style('font-size', '12px')
  .style('font-family', 'sans-serif')
  .style('fill', '#000000')
  .text(d => d.label);

  const present = (Object.keys(this.SEVERITY_STYLE) as Severity[])
  .filter(s => cols.some(c => c.severity === s));
if (present.length) {
  const sevLegend = mainGroup.append('g')
    .attr('transform', `translate(${legendX}, ${totalLegendHeight + 30})`);
  sevLegend.append('text').attr('y', -8)
    .style('font-size', '12px').style('font-weight', 'bold').text('Severity');
  present.forEach((s, i) => {
    const st = this.SEVERITY_STYLE[s];
    sevLegend.append('rect').attr('y', i * 18).attr('width', legendWidth).attr('height', 14)
      .attr('rx', 3).style('fill', st.fill).style('stroke', '#333333');
    sevLegend.append('text').attr('x', legendWidth + 10).attr('y', i * 18 + 7).attr('dy', '0.35em')
      .style('font-size', '12px').text(s.charAt(0).toUpperCase() + s.slice(1));
  });
  const noteY = present.length * 18 + 14;
  ['Thick outline =', 'matched severe term'].forEach((line, i) =>
    sevLegend.append('text').attr('y', noteY + i * 12)
      .style('font-size', '10px').style('fill', '#555555').text(line));
  const legendBottom = totalLegendHeight + 30 + noteY + 18;   // second note line plus descenders
  svg.attr('height', Math.max(height, margin.top + legendBottom));
}


  }

  
}


