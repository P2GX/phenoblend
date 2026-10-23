import {
  Component,
  ElementRef,
  input,
  Input,
  OnChanges,
  SimpleChanges,
  ViewChild,
} from '@angular/core';
import * as d3 from 'd3';
import { UpsetPlotPayload } from '../models/phenoprofile_dto';

@Component({
  selector: 'app-upset-plot',
  standalone: true,
  template: `<div #chartContainer class="upset-container"></div>`,
  styles: [
    `
      .upset-container {
        font-family:
          -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      }
    `,
  ],
})
export class UpsetPlotComponent implements OnChanges {
  @ViewChild('chartContainer', { static: true })
  chartContainerRef!: ElementRef<HTMLElement>;

  data = input.required<UpsetPlotPayload>();

  private readonly BASE_COLOR = '#8D10CB';
  private readonly HIGHLIGHT_COLOR = '#66C2A5';
  private readonly LIGHT_GREY = '#E0E0E0';

  ngOnChanges(changes: SimpleChanges): void {
    if (changes['data'] && this.data()) {
      this.renderUpsetPlot();
    }
  }

  private renderUpsetPlot(): void {
    const container = this.chartContainerRef.nativeElement;
    container.innerHTML = '';

    const {
      genes,
      combinations,
      combinationAnnotated,
      combinationObserved,
      geneAnnotated,
      geneObserved,
    } = this.data();

    const nGenes = genes.length;
    const nCombos = combinations.length;
    const cellSize = 35;
    const barHeightMax = 150;
    const marginalWidthMax = 150;

    const margin = { top: 30, right: 40, bottom: 40, left: 50 };

    const matrixWidth = nCombos * cellSize;
    const matrixHeight = nGenes * cellSize;
    const xMatrixOffset = margin.left + marginalWidthMax + 65; // 65px gap
    const yMatrixOffset = margin.top + barHeightMax + 10; // 10px gap
    // Complete view boundaries
    const totalWidth = xMatrixOffset + matrixWidth +  margin.right;
    const totalHeight = yMatrixOffset + matrixHeight + margin.bottom;

    const svg = d3
      .select(container)
      .append('svg')
      .attr('width', totalWidth)
      .attr('height', totalHeight);

    

    // ==========================================
    // 1. SCALES SETUP
    // ==========================================
    const xScale = d3
      .scaleBand()
      .domain(d3.range(nCombos).map(String))
      .range([0, matrixWidth])
      .padding(0.1);

    const yScale = d3
      .scaleBand()
      .domain(genes)
      .range([0, matrixHeight])
      .padding(0.1);

    const yBarScale = d3
      .scaleLinear()
      .domain([0, d3.max(combinationAnnotated) || 1])
      .range([barHeightMax, 0])
      .nice();

    const xMarginalScale = d3
      .scaleLinear()
      .domain([0, d3.max(geneAnnotated) || 1])
      .range([0, marginalWidthMax])
      .nice();

    // ==========================================
    // 2. TOP RIGHT: STACKED BAR CHART (ax_bar)
    // ==========================================
    const barGroup = svg
      .append('g')
      .attr('transform', `translate(${xMatrixOffset}, ${margin.top})`);

    // Annotated background bars
    barGroup
      .selectAll('.bar-annotated')
      .data(combinationAnnotated)
      .enter()
      .append('rect')
      .attr('x', (_, i) => xScale(String(i))!)
      .attr('y', (d) => yBarScale(d))
      .attr('width', xScale.bandwidth())
      .attr('height', (d) => barHeightMax - yBarScale(d))
      .style('fill', this.BASE_COLOR);

    // Observed foreground bars
    barGroup
      .selectAll('.bar-observed')
      .data(combinationObserved)
      .enter()
      .append('rect')
      .attr('x', (_, i) => xScale(String(i))!)
      .attr('y', (d) => yBarScale(d))
      .attr('width', xScale.bandwidth())
      .attr('height', (d) => barHeightMax - yBarScale(d))
      .style('fill', this.HIGHLIGHT_COLOR);

    // Y Axis for Intersections
    barGroup
      .append('g')
      .call(d3.axisLeft(yBarScale).ticks(5))
      .append('text')
      .attr('transform', 'rotate(-90)')
      .attr('y', -40)
      .attr('x', -barHeightMax / 2)
      .attr('fill', '#000')
      .style('text-anchor', 'middle')
      .style('font-size', '11px')
      .text('Disease models');

    // ==========================================
    // 3. BOTTOM LEFT: MARGINAL GENE BARS (ax_marginal)
    // ==========================================
    const marginalGroup = svg
      .append('g')
      .attr('transform', `translate(${margin.left}, ${yMatrixOffset})`);

    // Invert scale layout structure to grow leftwards (invert_xaxis)
    const reversedXMarginalScale = xMarginalScale
      .copy()
      .range([marginalWidthMax, 0]);

    // Total Annotated Bar lengths
    marginalGroup
      .selectAll('.marginal-annotated')
      .data(genes)
      .enter()
      .append('rect')
      .attr('x', (_, i) => reversedXMarginalScale(geneAnnotated[i]))
      .attr('y', (d) => yScale(d)!)
      .attr(
        'width',
        (_, i) => marginalWidthMax - reversedXMarginalScale(geneAnnotated[i]),
      )
      .attr('height', yScale.bandwidth())
      .style('fill', this.BASE_COLOR);

    // Overlapping Observed Bar lengths
    marginalGroup
      .selectAll('.marginal-observed')
      .data(genes)
      .enter()
      .append('rect')
      .attr('x', (_, i) => reversedXMarginalScale(geneObserved[i]))
      .attr('y', (d) => yScale(d)!)
      .attr(
        'width',
        (_, i) => marginalWidthMax - reversedXMarginalScale(geneObserved[i]),
      )
      .attr('height', yScale.bandwidth())
      .style('fill', this.HIGHLIGHT_COLOR);

    // X Axis at the bottom of the marginal plot
    marginalGroup
      .append('g')
      .attr('transform', `translate(0, ${matrixHeight})`)
      .call(d3.axisBottom(reversedXMarginalScale).ticks(4))
      .append('text')
      .attr('x', marginalWidthMax / 2)
      .attr('y', 30)
      .attr('fill', '#000')
      .style('text-anchor', 'middle')
      .style('font-size', '11px')
      .text('Individual');

    // ==========================================
    // 4. BOTTOM RIGHT: DOT MATRIX LAYOUT (ax_matrix)
    // ==========================================
    const matrixGroup = svg
      .append('g')
      .attr('transform', `translate(${xMatrixOffset}, ${yMatrixOffset})`);

    // Render underlying layout connection spans
    combinations.forEach((combo, colIdx) => {
      const activeRows = combo.map(
        (gene) => yScale(gene)! + yScale.bandwidth() / 2,
      );
      if (activeRows.length > 1) {
        matrixGroup
          .append('line')
          .attr('x1', xScale(String(colIdx))! + xScale.bandwidth() / 2)
          .attr('x2', xScale(String(colIdx))! + xScale.bandwidth() / 2)
          .attr('y1', d3.min(activeRows)!)
          .attr('y2', d3.max(activeRows)!)
          .style('stroke', this.BASE_COLOR)
          .style('stroke-width', '2.5px');
      }
    });

    // Populate discrete nodes grid layout
    combinations.forEach((combo, colIdx) => {
      genes.forEach((gene) => {
        const isPresent = combo.includes(gene);
        matrixGroup
          .append('circle')
          .attr('cx', xScale(String(colIdx))! + xScale.bandwidth() / 2)
          .attr('cy', yScale(gene)! + yScale.bandwidth() / 2)
          .attr('r', 6)
          .style('fill', isPresent ? this.BASE_COLOR : this.LIGHT_GREY);
      });
    });

    // Y Axis Labels for individual Genes
    matrixGroup
      .append('g')
      .call(d3.axisLeft(yScale).tickSize(0))
      .select('.domain')
      .remove();

    matrixGroup
      .selectAll('.tick text')
      .style('font-size', '12px')
      .style('font-weight', 'bold');

    // ==========================================
    // 5. TOP LEFT: LEGEND AREA (ax_legend)
    // ==========================================
    const legendGroup = svg
      .append('g')
      .attr('transform', `translate(${margin.left}, ${margin.top + 20})`);

    const legendItems = [
      { text: 'Annotated', color: this.BASE_COLOR },
      { text: 'Observed', color: this.HIGHLIGHT_COLOR },
    ];

    legendItems.forEach((item, idx) => {
      const g = legendGroup
        .append('g')
        .attr('transform', `translate(0, ${idx * 20})`);

      g.append('rect')
        .attr('width', 15)
        .attr('height', 15)
        .style('fill', item.color);

      g.append('text')
        .attr('x', 22)
        .attr('y', 12)
        .style('font-size', '12px')
        .text(item.text);
    });
  }
}

export { UpsetPlotPayload };
