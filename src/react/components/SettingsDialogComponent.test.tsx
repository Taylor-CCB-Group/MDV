import { describe, test, expect, vi } from 'vitest';
import { g } from '@/lib/utils';
import type { AnyGuiSpec } from '@/charts/charts';
import { collectDisposers } from './SettingsDialogComponent';
import type { IReactionDisposer } from 'mobx';

// We need to export collectDisposers for testing
// Since it's not exported, we'll test it indirectly through the component
// or we can make it exported. For now, let's test the disposer collection logic directly.

describe('Settings Dialog Disposer Management', () => {
    test('collectDisposers collects disposers from flat specs', () => {
        const disposer1 = vi.fn() as unknown as IReactionDisposer;
        const disposer2 = vi.fn() as unknown as IReactionDisposer;
        
        const specs: AnyGuiSpec[] = [
            g({
                type: 'slider',
                label: 'Test Slider 1',
                current_value: 0.5,
                _disposers: [disposer1],
            }),
            g({
                type: 'slider',
                label: 'Test Slider 2',
                current_value: 0.7,
                _disposers: [disposer2],
            }),
        ];
        
        const collected = collectDisposers(specs);
        expect(collected).toHaveLength(2);
        expect(collected).toContain(disposer1);
        expect(collected).toContain(disposer2);
    });
    
    test('collectDisposers collects disposers from nested folders', () => {
        const disposer1 = vi.fn() as unknown as IReactionDisposer;
        const disposer2 = vi.fn() as unknown as IReactionDisposer;
        const disposer3 = vi.fn() as unknown as IReactionDisposer;
        
        const specs: AnyGuiSpec[] = [
            g({
                type: 'folder',
                label: 'Outer Folder',
                current_value: [
                    g({
                        type: 'slider',
                        label: 'Nested Slider',
                        current_value: 0.5,
                        _disposers: [disposer1],
                    }),
                    g({
                        type: 'folder',
                        label: 'Inner Folder',
                        current_value: [
                            g({
                                type: 'slider',
                                label: 'Deep Slider',
                                current_value: 0.3,
                                _disposers: [disposer2, disposer3],
                            }),
                        ],
                    }),
                ],
            }),
        ];
        
        const collected = collectDisposers(specs);
        expect(collected).toHaveLength(3);
        expect(collected).toContain(disposer1);
        expect(collected).toContain(disposer2);
        expect(collected).toContain(disposer3);
    });
    
    test('collectDisposers handles specs without disposers', () => {
        const disposer1 = vi.fn() as unknown as IReactionDisposer;
        
        const specs: AnyGuiSpec[] = [
            g({
                type: 'slider',
                label: 'Slider with disposer',
                current_value: 0.5,
                _disposers: [disposer1],
            }),
            g({
                type: 'slider',
                label: 'Slider without disposer',
                current_value: 0.7,
            }),
        ];
        
        const collected = collectDisposers(specs);
        expect(collected).toHaveLength(1);
        expect(collected).toContain(disposer1);
    });
    
    test('collectDisposers handles empty specs array', () => {
        const specs: AnyGuiSpec[] = [];
        const collected = collectDisposers(specs);
        expect(collected).toHaveLength(0);
    });
    
    test('collectDisposers handles empty folders', () => {
        const specs: AnyGuiSpec[] = [
            g({
                type: 'folder',
                label: 'Empty Folder',
                current_value: [],
            }),
        ];
        
        const collected = collectDisposers(specs);
        expect(collected).toHaveLength(0);
    });
    
    test('collectDisposers collects multiple disposers from same spec', () => {
        const disposer1 = vi.fn() as unknown as IReactionDisposer;
        const disposer2 = vi.fn() as unknown as IReactionDisposer;
        const disposer3 = vi.fn() as unknown as IReactionDisposer;
        
        const specs: AnyGuiSpec[] = [
            g({
                type: 'slider',
                label: 'Slider with multiple disposers',
                current_value: 0.5,
                _disposers: [disposer1, disposer2, disposer3],
            }),
        ];
        
        const collected = collectDisposers(specs);
        expect(collected).toHaveLength(3);
        expect(collected).toContain(disposer1);
        expect(collected).toContain(disposer2);
        expect(collected).toContain(disposer3);
    });
    
    test('collectDisposers handles deeply nested folder structures', () => {
        const disposers: IReactionDisposer[] = [];
        for (let i = 0; i < 5; i++) {
            disposers.push(vi.fn() as unknown as IReactionDisposer);
        }
        
        const specs: AnyGuiSpec[] = [
            g({
                type: 'folder',
                label: 'Level 1',
                current_value: [
                    g({
                        type: 'folder',
                        label: 'Level 2',
                        current_value: [
                            g({
                                type: 'folder',
                                label: 'Level 3',
                                current_value: [
                                    g({
                                        type: 'slider',
                                        label: 'Deep Slider',
                                        current_value: 0.5,
                                        _disposers: [disposers[0], disposers[1]],
                                    }),
                                ],
                                _disposers: [disposers[2]],
                            }),
                        ],
                        _disposers: [disposers[3]],
                    }),
                ],
                _disposers: [disposers[4]],
            }),
        ];
        
        const collected = collectDisposers(specs);
        expect(collected).toHaveLength(5);
        disposers.forEach((disposer) => {
            expect(collected).toContain(disposer);
        });
    });
});

