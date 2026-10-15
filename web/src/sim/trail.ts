/** Fixed-capacity ring buffer of xyz points (oldest overwritten). */
export class RingTrail {
  readonly capacity: number;
  private buf: Float32Array;
  private head = 0; // next write index
  count = 0;
  constructor(capacity: number) {
    this.capacity = capacity;
    this.buf = new Float32Array(3 * capacity);
  }
  clear(): void { this.head = 0; this.count = 0; }
  push(x: number, y: number, z: number): void {
    const i = 3 * this.head;
    this.buf[i] = x; this.buf[i + 1] = y; this.buf[i + 2] = z;
    this.head = (this.head + 1) % this.capacity;
    if (this.count < this.capacity) this.count++;
  }
  /** Copy points oldest-first into out; returns count. */
  copyOrdered(out: Float32Array): number {
    const start = (this.head - this.count + this.capacity) % this.capacity;
    for (let k = 0; k < this.count; k++) {
      const s = 3 * ((start + k) % this.capacity);
      out[3 * k] = this.buf[s]!; out[3 * k + 1] = this.buf[s + 1]!; out[3 * k + 2] = this.buf[s + 2]!;
    }
    return this.count;
  }
}
