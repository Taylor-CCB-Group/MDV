/**
 * MDV's data loaders need SharedArrayBuffer, which a browser provides only to a cross-origin
 * isolated page. A page inside a frame is isolated only when the page around it is too, so a
 * host that frames MDV without isolation, such as Verily Workbench, breaks data loading. In a
 * tab of its own the same MDV works, so in that case MDV offers one.
 */

export type FrameState = { top: unknown; crossOriginIsolated: boolean };

const NOTICE_ID = "mdv-own-tab-notice";

export function needsOwnTab(win: FrameState = window): boolean {
    return win.top !== win && !win.crossOriginIsolated;
}

export function showOwnTabNotice(doc: Document, href: string): void {
    if (doc.getElementById(NOTICE_ID)) return;

    const notice = doc.createElement("div");
    notice.id = NOTICE_ID;
    notice.setAttribute("role", "alert");
    notice.style.cssText =
        "position:fixed;top:0;left:0;right:0;z-index:2147483647;padding:10px 16px;" +
        "background:#fff4e5;color:#663c00;border-bottom:1px solid #ffb74d;" +
        "font:14px/1.4 system-ui,sans-serif;";

    const link = doc.createElement("a");
    link.href = href;
    link.target = "_blank";
    link.rel = "noopener";
    link.textContent = "Open MDV in a new tab";
    link.style.fontWeight = "600";

    notice.append(
        "MDV cannot load data inside this frame, because the page around it does not allow it. ",
        link,
        ` to use it, or open ${href} directly.`,
    );
    doc.body.prepend(notice);
}

export function offerOwnTabIfFramed(): void {
    if (needsOwnTab()) showOwnTabNotice(document, window.location.href);
}
