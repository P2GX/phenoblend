// notebook.component.ts
import { Component, computed, inject } from '@angular/core';
import { DatePipe } from '@angular/common';
import { AppStatusService } from '../services/app-status-service';

@Component({
  selector: 'app-notebook',
  standalone: true,
  imports: [DatePipe],
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
}