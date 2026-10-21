import { Meta, StoryObj, moduleMetadata } from '@storybook/angular';
import { CommonModule, DecimalPipe } from '@angular/common';
import { OverlapPlotComponent } from './overlap-plot.component';

const meta: Meta<OverlapPlotComponent> = {
  title: 'Phenoblend/OverlapPlotPayload',
  component: OverlapPlotComponent,
  decorators: [
    moduleMetadata({
      imports: [CommonModule, DecimalPipe],
    }),
  ],
  tags: ['autodocs'],
};

export default meta;
type Story = StoryObj<OverlapPlotComponent>;

// Sample mock payload mimicking your D3 payload DTO
const mockPayload = {
  entities: ['FGFR2', 'TWIST1', 'SOX9'],
  columns: [
    {
      hpoId: 'HP:0000248',
      hpoName: 'Brachycephaly',
      scores: { 'FGFR2': 1.0, 'TWIST1': 0.8, 'SOX9': 0.0 }
    },
    {
      hpoId: 'HP:0000218',
      hpoName: 'High forehead',
      scores: { 'FGFR2': 0.5, 'TWIST1': 1.0, 'SOX9': 0.2 }
    },
    {
      hpoId: 'HP:0001249',
      hpoName: 'Intellectual disability but with long wrap text testing',
      scores: { 'FGFR2': 0.0, 'TWIST1': 0.1, 'SOX9': 0.9 }
    }
  ]
};

export const Default: Story = {
  args: {
    data: mockPayload,
  },
};

export const LargeGenePanel: Story = {
  args: {
    data: {
      entities: Array.from({ length: 15 }, (_, i) => `GENE-${i + 1}`),
      columns: mockPayload.columns,
    },
  },
};  