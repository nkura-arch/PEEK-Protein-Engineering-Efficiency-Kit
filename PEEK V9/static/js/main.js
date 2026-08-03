// Main JavaScript for PEEK application

/* ============================================================
 * Cross-page input persistence
 * ------------------------------------------------------------
 * Saves form field values to localStorage so the DNA sequence
 * and other inputs are NOT lost when navigating between pages
 * (Silent Mutations <-> Primer Design <-> Sequence Visualizer).
 * Fields are keyed by their `name`, so the shared "dna" field
 * carries across every tab automatically.
 * ============================================================ */
const PEEK_STORE_PREFIX = 'peek:field:';
// Field names that should NOT be persisted (transient / per-submit only).
const PEEK_NO_PERSIST = new Set(['csrf_token']);

function peekPersistKey(name) {
    return PEEK_STORE_PREFIX + name;
}

function peekSaveField(el) {
    if (!el.name || PEEK_NO_PERSIST.has(el.name)) return;
    try {
        if (el.type === 'checkbox') {
            localStorage.setItem(peekPersistKey(el.name), el.checked ? '1' : '0');
        } else if (el.type === 'radio') {
            // Only the selected radio in a group should record the group's value,
            // otherwise the last radio in the DOM always wins.
            if (el.checked) localStorage.setItem(peekPersistKey(el.name), el.value);
        } else {
            localStorage.setItem(peekPersistKey(el.name), el.value);
        }
    } catch (e) { /* storage full / disabled — ignore */ }
}

function peekRestoreFields() {
    // On a plain navigation (GET) the server only renders defaults, so a saved
    // value must win even over a non-empty default (e.g. Max Changes = 3) —
    // otherwise every field with a default silently resets between tabs.
    // On a POST the server is echoing the values just submitted, so we only
    // fill fields that are genuinely empty and never clobber the real result.
    const navRestore = !window.PEEK_DID_POST;
    const controls = document.querySelectorAll('input[name], textarea[name], select[name]');
    controls.forEach(el => {
        if (!el.name || PEEK_NO_PERSIST.has(el.name)) return;
        if (el.type === 'submit' || el.type === 'button' || el.type === 'file' ||
            el.type === 'hidden' || el.type === 'password') return;

        let stored;
        try { stored = localStorage.getItem(peekPersistKey(el.name)); }
        catch (e) { return; }
        if (stored === null) return;

        if (el.type === 'checkbox') {
            if (navRestore) el.checked = stored === '1';
        } else if (el.type === 'radio') {
            // Restore which option in the group was selected.
            if (navRestore && el.value === stored && !el.checked) {
                el.checked = true;
                el.dispatchEvent(new Event('change', { bubbles: true }));
            }
        } else if (el.tagName === 'SELECT') {
            if ((navRestore || !el.value || el.selectedIndex < 0) && el.value !== stored) {
                el.value = stored;
                el.dispatchEvent(new Event('change', { bubbles: true }));
            }
        } else {
            const shouldSet = navRestore
                ? (stored.trim() !== '' && el.value !== stored)
                : (el.value.trim() === '' && stored.trim() !== '');
            if (shouldSet) {
                el.value = stored;
                // Notify listeners (bp counters, validation, uppercase, etc.)
                el.dispatchEvent(new Event('input', { bubbles: true }));
            }
        }
    });
}

function peekInitPersistence() {
    peekRestoreFields();
    // After restoring saved input, offer a default example on untouched fields.
    peekApplyExample();
    // Save on edit.
    document.addEventListener('input', function(e) {
        const t = e.target;
        if (t && (t.matches('input[name], textarea[name], select[name]'))) {
            peekSaveField(t);
        }
    });
    document.addEventListener('change', function(e) {
        const t = e.target;
        if (t && (t.matches('input[name], textarea[name], select[name]'))) {
            peekSaveField(t);
        }
    });
    // Save everything right before submitting (captures programmatic fills).
    document.addEventListener('submit', function(e) {
        if (e.target && e.target.tagName === 'FORM') {
            e.target.querySelectorAll('input[name], textarea[name], select[name]')
                .forEach(peekSaveField);
        }
    }, true);
}

// Clear all persisted PEEK fields (exposed for an optional "reset" button).
function peekClearSaved() {
    try {
        Object.keys(localStorage)
            .filter(k => k.startsWith(PEEK_STORE_PREFIX) || k.indexOf('peek:ran:') === 0)
            .forEach(k => localStorage.removeItem(k));
    } catch (e) { /* ignore */ }
    try {
        Object.keys(sessionStorage)
            .filter(k => k.indexOf('peek:res:') === 0 || k.indexOf('peek:resreplay:') === 0)
            .forEach(k => sessionStorage.removeItem(k));
    } catch (e) { /* ignore */ }
}

/* ============================================================
 * Default example loader
 * ------------------------------------------------------------
 * A page may set `window.PEEK_EXAMPLE = { fieldName: value, ... }`.
 * peekApplyExample() fills a field ONLY when the user has never
 * touched it (no persisted value) and it is currently empty, so
 * real input and server-rendered results are never overwritten.
 * The shared "dna" field uses the same example across tabs, so an
 * example sequence carries from one tab to the next just like
 * user input does.
 * ============================================================ */
function peekFillField(el, value) {
    if (!el) return;
    el.value = value;
    // Notify listeners (bp counters, uppercasing, validation, persistence).
    el.dispatchEvent(new Event('input', { bubbles: true }));
    el.dispatchEvent(new Event('change', { bubbles: true }));
}

function peekApplyExample() {
    const ex = window.PEEK_EXAMPLE;
    if (!ex) return;
    Object.keys(ex).forEach(function (name) {
        if (name.charAt(0) === '_') return;
        const el = document.querySelector('[name="' + name + '"]');
        if (!el) return;
        let stored;
        try { stored = localStorage.getItem(peekPersistKey(name)); } catch (e) { return; }
        if (stored !== null) return;                  // user/app already set this field
        if ((el.value || '').trim() !== '') return;   // server already rendered a value
        peekFillField(el, ex[name]);
    });
}

// Force-load the page example (used by a "Load example" button).
function peekLoadExample() {
    const ex = window.PEEK_EXAMPLE;
    if (!ex) return;
    Object.keys(ex).forEach(function (name) {
        if (name.charAt(0) === '_') return;
        const el = document.querySelector('[name="' + name + '"]');
        if (el) peekFillField(el, ex[name]);
    });
}

// Clear every named field on the page and forget the saved values.
function peekClearForm() {
    const controls = document.querySelectorAll('input[name], textarea[name], select[name]');
    controls.forEach(function (el) {
        if (!el.name || PEEK_NO_PERSIST.has(el.name)) return;
        if (el.type === 'submit' || el.type === 'button' || el.type === 'hidden') return;
        if (el.type === 'checkbox') {
            // leave checkboxes at their current (default) state
        } else if (el.tagName === 'SELECT') {
            el.selectedIndex = 0;
            el.dispatchEvent(new Event('change', { bubbles: true }));
        } else {
            el.value = '';
            el.dispatchEvent(new Event('input', { bubbles: true }));
        }
        try { localStorage.removeItem(peekPersistKey(el.name)); } catch (e) { /* ignore */ }
    });
    // Drop cached results for a page the user just cleared.
    const key = peekCurrentPageKey();
    if (key) peekClearResultsCache(key);
}

document.addEventListener('DOMContentLoaded', peekInitPersistence);

/* ============================================================
 * Results persistence across navigation
 * ------------------------------------------------------------
 * Goal: a page's results STAY when you switch to another subtab,
 * and reset only when you (a) run another iteration, or (b) refresh
 * the tab. Inputs are kept either way.
 *
 * Server-rendered results only appear on POST, so a plain GET back
 * to a tab would normally lose them. We solve this WITHOUT re-running:
 * on a successful POST we cache the rendered results HTML in
 * sessionStorage (per tab) and restore it instantly on return.
 * sessionStorage is per-tab, so results never leak between browser
 * tabs and vanish when the tab is closed.
 *
 * Refresh detection: the Navigation Timing API tells us whether this
 * load was a reload vs a normal navigation. On a reload we drop the
 * cached results so a refresh starts clean (inputs still persist).
 *
 * Fallback: a single run can be thousands of rows / many MB, which
 * can exceed the storage quota. When caching fails we fall back to
 * replaying the (already-restored) input form on return.
 *
 * A page opts in by wrapping its results-only markup in:
 *     <div id="peekResults" data-peek-page="KEY"> ... </div>
 * After results are present, the page may define
 * window.peekOnResultsReady(key) to re-init scripted widgets
 * (e.g. an SVG drawn from a data island).
 * ============================================================ */
function peekResHtmlKey(key)  { return 'peek:res:' + key; }      // cached results HTML
function peekResReplayKey(key) { return 'peek:resreplay:' + key; } // "too big — replay" marker

function peekAfterResults(key) {
    if (typeof window.peekOnResultsReady === 'function') {
        try { window.peekOnResultsReady(key); } catch (e) { /* ignore */ }
    }
}

// The primary input form is the first <form> inside <main> on every tab.
function peekMainForm() {
    return document.querySelector('main form') || document.querySelector('form');
}

function peekCurrentPageKey() {
    const wrap = document.getElementById('peekResults');
    return wrap ? (wrap.getAttribute('data-peek-page') || 'default') : null;
}

// 'reload' | 'navigate' | 'back_forward' | 'navigate' (fallback).
function peekNavType() {
    try {
        const e = (performance.getEntriesByType('navigation') || [])[0];
        if (e && e.type) return e.type;
    } catch (e) { /* ignore */ }
    try {
        if (performance.navigation) {
            const t = performance.navigation.type;
            if (t === 1) return 'reload';
            if (t === 2) return 'back_forward';
        }
    } catch (e) { /* ignore */ }
    return 'navigate';
}

function peekClearResultsCache(key) {
    try { sessionStorage.removeItem(peekResHtmlKey(key)); } catch (e) { /* ignore */ }
    try { sessionStorage.removeItem(peekResReplayKey(key)); } catch (e) { /* ignore */ }
}

function peekInitResultsPersistence() {
    // Clean up keys written by previous versions of this code.
    try {
        Object.keys(localStorage)
            .filter(k => k.indexOf('peek:results:') === 0 || k.indexOf('peek:ran:') === 0)
            .forEach(k => localStorage.removeItem(k));
    } catch (e) { /* ignore */ }

    const wrap = document.getElementById('peekResults');
    if (!wrap) return;
    const key = wrap.getAttribute('data-peek-page') || 'default';
    const emptyState = document.getElementById('peekEmptyState');
    // The server rendered results into the wrapper iff it has element children.
    const hasServerResults = wrap.children.length > 0;
    const isReload = peekNavType() === 'reload';

    // (a) Refresh resets results: forget the cache for this tab. If the reload
    //     itself rendered results (e.g. a resubmitted POST), keep showing them.
    if (isReload) {
        peekClearResultsCache(key);
        if (hasServerResults) peekAfterResults(key);
        return;
    }

    // (b) A run just happened (POST): the server output is the new truth.
    if (window.PEEK_DID_POST) {
        if (hasServerResults) {
            let cached = false;
            try { sessionStorage.setItem(peekResHtmlKey(key), wrap.innerHTML); cached = true; }
            catch (e) { cached = false; }
            try { sessionStorage.removeItem(peekResReplayKey(key)); } catch (e) { /* ignore */ }
            // Too large to cache — remember to replay the form on return instead.
            if (!cached) {
                try { sessionStorage.setItem(peekResReplayKey(key), '1'); } catch (e) { /* ignore */ }
            }
            peekAfterResults(key);
        } else {
            peekClearResultsCache(key);
        }
        return;
    }

    // (c) Plain navigation back to this tab: restore what was there.
    if (hasServerResults) { peekAfterResults(key); return; }

    let cachedHtml = null;
    try { cachedHtml = sessionStorage.getItem(peekResHtmlKey(key)); } catch (e) { cachedHtml = null; }
    if (cachedHtml) {
        wrap.innerHTML = cachedHtml;
        if (emptyState) emptyState.style.display = 'none';
        peekAfterResults(key);
        return;
    }

    let needReplay = null;
    try { needReplay = sessionStorage.getItem(peekResReplayKey(key)); } catch (e) { needReplay = null; }
    if (needReplay !== '1') return;

    const form = peekMainForm();
    if (!form) return;
    // Only replay once the required inputs are present (restored above);
    // otherwise leave the empty state so we never loop on a blank form.
    let ready = true;
    form.querySelectorAll('[required]').forEach(el => {
        if (!String(el.value || '').trim()) ready = false;
    });
    if (!ready) return;
    if (emptyState) emptyState.style.display = 'none'; // reduce flash before re-render
    if (typeof form.requestSubmit === 'function') form.requestSubmit();
    else form.submit();
}

document.addEventListener('DOMContentLoaded', peekInitResultsPersistence);

// Copy to clipboard functionality
document.addEventListener('DOMContentLoaded', function() {
    // Add click handlers to copy buttons
    document.addEventListener('click', function(e) {
        if (e.target.classList.contains('copy-btn') || e.target.closest('.copy-btn')) {
            const button = e.target.classList.contains('copy-btn') ? e.target : e.target.closest('.copy-btn');
            const textToCopy = button.getAttribute('data-copy');
            
            if (textToCopy) {
                copyToClipboard(textToCopy);
                showCopySuccess(button);
            }
        }
    });

    // Form validation
    const forms = document.querySelectorAll('form');
    forms.forEach(form => {
        form.addEventListener('submit', function(e) {
            console.log('Form submit attempt:', form.action, form.method);
            const isValid = validateForm(form);
            console.log('Form validation result:', isValid);
            
            if (!isValid) {
                console.log('Form validation failed, preventing submission');
                e.preventDefault();
                return false;
            }
            
            console.log('Form is valid, proceeding with submission');
            // Add loading state
            const submitBtn = form.querySelector('button[type="submit"]');
            if (submitBtn) {
                addLoadingState(submitBtn);
            }
        });
    });

    // DNA sequence validation
    const dnaInputs = document.querySelectorAll('textarea[name="dna"], input[name="dna"]');
    dnaInputs.forEach(input => {
        input.addEventListener('input', function(e) {
            validateDNASequence(e.target);
        });
    });

    // Auto-uppercase DNA sequences
    const sequenceInputs = document.querySelectorAll('textarea[name="dna"], textarea[name="aa_seq"]');
    sequenceInputs.forEach(input => {
        input.addEventListener('input', function(e) {
            e.target.value = e.target.value.toUpperCase();
        });
    });
});

/**
 * Copy text to clipboard
 */
async function copyToClipboard(text) {
    try {
        await navigator.clipboard.writeText(text);
        return true;
    } catch (err) {
        // Fallback for older browsers
        const textArea = document.createElement('textarea');
        textArea.value = text;
        document.body.appendChild(textArea);
        textArea.focus();
        textArea.select();
        try {
            document.execCommand('copy');
            document.body.removeChild(textArea);
            return true;
        } catch (fallbackErr) {
            document.body.removeChild(textArea);
            return false;
        }
    }
}

/**
 * Show copy success feedback
 */
function showCopySuccess(button) {
    const originalContent = button.innerHTML;
    button.innerHTML = '<i class="fas fa-check"></i> Copied!';
    button.classList.add('btn-success');
    button.classList.remove('btn-outline-primary', 'btn-outline-light');
    
    setTimeout(() => {
        button.innerHTML = originalContent;
        button.classList.remove('btn-success');
        button.classList.add('btn-outline-primary');
    }, 2000);
}

/**
 * Validate DNA sequence
 */
function validateDNASequence(input) {
    const sequence = input.value.toUpperCase();
    const validBases = /^[ATCG]*$/;
    
    if (sequence && !validBases.test(sequence)) {
        input.setCustomValidity('DNA sequence must contain only A, T, C, G characters');
        input.classList.add('is-invalid');
        
        // Add or update error message
        let errorDiv = input.parentNode.querySelector('.invalid-feedback');
        if (!errorDiv) {
            errorDiv = document.createElement('div');
            errorDiv.className = 'invalid-feedback';
            input.parentNode.appendChild(errorDiv);
        }
        errorDiv.textContent = 'DNA sequence must contain only A, T, C, G characters';
    } else {
        input.setCustomValidity('');
        input.classList.remove('is-invalid');
        
        // Remove error message
        const errorDiv = input.parentNode.querySelector('.invalid-feedback');
        if (errorDiv) {
            errorDiv.remove();
        }
    }
}

/**
 * Validate form before submission
 */
function validateForm(form) {
    let isValid = true;
    console.log('Validating form...');
    
    // Validate required fields
    const requiredFields = form.querySelectorAll('[required]');
    console.log('Required fields:', requiredFields.length);
    requiredFields.forEach(field => {
        if (!field.value.trim()) {
            console.log('Required field empty:', field.name);
            field.classList.add('is-invalid');
            isValid = false;
        } else {
            console.log('Required field valid:', field.name, field.value.length, 'chars');
            field.classList.remove('is-invalid');
        }
    });
    
    // Validate DNA sequence
    const dnaField = form.querySelector('textarea[name="dna"]');
    if (dnaField && dnaField.value) {
        const sequence = dnaField.value.toUpperCase().replace(/\s+/g, ''); // Remove whitespace
        const validBases = /^[ATCG]*$/;
        console.log('DNA sequence length:', sequence.length, 'valid:', validBases.test(sequence));
        if (!validBases.test(sequence)) {
            dnaField.classList.add('is-invalid');
            console.log('DNA validation failed');
            isValid = false;
        } else {
            dnaField.classList.remove('is-invalid');
        }
    }
    
    // Validate amino acid sequence if present
    const aaField = form.querySelector('textarea[name="aa_seq"]');
    if (aaField && aaField.value) {
        const aaSequence = aaField.value.toUpperCase();
        const validAA = /^[ARNDCEQGHILKMFPSTWYV*]*$/;
        if (!validAA.test(aaSequence)) {
            aaField.classList.add('is-invalid');
            showError('Amino acid sequence contains invalid characters');
            isValid = false;
        }
    }
    
    console.log('Form validation complete, isValid:', isValid);
    return isValid;
}

/**
 * Add loading state to button
 */
function addLoadingState(button) {
    const originalText = button.innerHTML;
    button.innerHTML = '<span class="spinner"></span> Processing...';
    button.disabled = true;
    
    // Reset after 30 seconds max
    setTimeout(() => {
        button.innerHTML = originalText;
        button.disabled = false;
    }, 30000);
}

/**
 * Show error message
 */
function showError(message) {
    const alertDiv = document.createElement('div');
    alertDiv.className = 'alert alert-danger alert-dismissible fade show';
    alertDiv.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    `;
    
    const container = document.querySelector('.container');
    if (container) {
        container.insertBefore(alertDiv, container.firstChild);
        
        // Auto-dismiss after 5 seconds
        setTimeout(() => {
            if (alertDiv.parentNode) {
                alertDiv.remove();
            }
        }, 5000);
    }
}

/**
 * Show success message
 */
function showSuccess(message) {
    const alertDiv = document.createElement('div');
    alertDiv.className = 'alert alert-success alert-dismissible fade show';
    alertDiv.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    `;
    
    const container = document.querySelector('.container');
    if (container) {
        container.insertBefore(alertDiv, container.firstChild);
        
        // Auto-dismiss after 3 seconds
        setTimeout(() => {
            if (alertDiv.parentNode) {
                alertDiv.remove();
            }
        }, 3000);
    }
}

/**
 * Format sequence with spaces for better readability
 */
function formatSequence(sequence, spacer = 10) {
    return sequence.replace(new RegExp(`(.{${spacer}})`, 'g'), '$1 ').trim();
}

/**
 * Calculate basic sequence statistics
 */
function calculateSequenceStats(sequence) {
    const length = sequence.length;
    const gcCount = (sequence.match(/[GC]/g) || []).length;
    const gcContent = ((gcCount / length) * 100).toFixed(1);
    
    return {
        length: length,
        gcContent: gcContent + '%',
        atContent: (100 - gcContent).toFixed(1) + '%'
    };
}

// Export functions for use in templates
window.PEEK = {
    copyToClipboard,
    showCopySuccess,
    validateDNASequence,
    formatSequence,
    calculateSequenceStats,
    showError,
    showSuccess,
    loadExample: peekLoadExample,
    clearForm: peekClearForm
};
