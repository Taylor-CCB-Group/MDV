const lazyChartRegistrationLoaders: Record<string, () => Promise<unknown>> = {
    table_chart_react: () => import("../react/components/TableChartReactWrapper"),
};

const lazyChartRegistrationPromises = new Map<string, Promise<unknown>>();

export function ensureLazyChartTypeRegistered(type: string): Promise<unknown> | null {
    const loader = lazyChartRegistrationLoaders[type];
    if (!loader) {
        return null;
    }
    const existing = lazyChartRegistrationPromises.get(type);
    if (existing) {
        return existing;
    }
    const promise = loader();
    lazyChartRegistrationPromises.set(type, promise);
    return promise;
}

export async function preloadLazyChartTypesForPicker() {
    await Promise.all(
        Object.keys(lazyChartRegistrationLoaders).map((type) =>
            ensureLazyChartTypeRegistered(type),
        ),
    );
}
