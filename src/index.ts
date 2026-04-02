import { shouldRenderDashboard } from "@/utils/mdvRouting";

async function bootstrap() {
    if (shouldRenderDashboard()) {
        document.getElementById("holder")?.remove();
        await import("./catalog/catalog_index");
        return;
    }

    await import("./modules/static_index");
}

void bootstrap();
