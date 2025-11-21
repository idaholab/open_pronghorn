(() => {
  const onReady = (cb) => {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', cb, { once: true });
    } else {
      cb();
    }
  };

  const parseDurationSeconds = (value, fallback = 18) => {
    if (!value) {
      return fallback;
    }
    const trimmed = value.toString().trim();
    if (!trimmed) {
      return fallback;
    }
    if (trimmed.endsWith('ms')) {
      return parseFloat(trimmed) / 1000;
    }
    return parseFloat(trimmed.replace('s', ''));
  };

  onReady(() => {
    const belt = document.querySelector('.logo-belt');
    if (!belt) {
      return;
    }

    const track = belt.querySelector('.logo-belt__track');
    if (!track || track.children.length < 2) {
      return;
    }

    const prefersReduced = window.matchMedia('(prefers-reduced-motion: reduce)');
    if (prefersReduced.matches) {
      return;
    }

    const items = Array.from(track.children);
    const state = {
      gap: 0,
      pxPerMs: 0,
      offset: 0,
      lastTs: null,
      rafId: null,
      paused: false,
    };

    const computeGap = () => {
      const style = getComputedStyle(track);
      const columnGap = parseFloat(style.columnGap || style.gap || '0');
      return Number.isNaN(columnGap) ? 0 : columnGap;
    };

    const updateMetrics = () => {
      state.gap = computeGap();
      const customDuration =
        belt.style.getPropertyValue('--logo-speed') ||
        getComputedStyle(belt).getPropertyValue('--logo-speed');
      const durationSeconds = parseDurationSeconds(customDuration) || 18;
      const totalWidth = items.reduce((sum, node, index) => {
        const rect = node.getBoundingClientRect();
        const gapContribution = index === items.length - 1 ? 0 : state.gap;
        return sum + rect.width + gapContribution;
      }, 0);
      state.pxPerMs = durationSeconds > 0 ? totalWidth / (durationSeconds * 1000) : 0;
    };

    const widthOfFirst = () => {
      const first = track.firstElementChild;
      if (!first) {
        return 0;
      }
      const rect = first.getBoundingClientRect();
      return rect.width + state.gap;
    };

    const applyTransform = () => {
      track.style.transform = `translateX(${-state.offset}px)`;
    };

    const advance = (deltaMs) => {
      if (state.pxPerMs <= 0) {
        return;
      }
      state.offset += state.pxPerMs * deltaMs;
      let firstWidth = widthOfFirst();
      while (firstWidth > 0 && state.offset >= firstWidth) {
        state.offset -= firstWidth;
        track.appendChild(track.firstElementChild);
        firstWidth = widthOfFirst();
      }
      applyTransform();
    };

    const step = (timestamp) => {
      if (state.lastTs === null) {
        state.lastTs = timestamp;
      }
      if (!state.paused) {
        const delta = timestamp - state.lastTs;
        advance(delta);
      }
      state.lastTs = timestamp;
      state.rafId = requestAnimationFrame(step);
    };

    const startLoop = () => {
      if (state.rafId !== null) {
        return;
      }
      state.lastTs = null;
      state.rafId = requestAnimationFrame(step);
    };

    const stopLoop = () => {
      if (state.rafId !== null) {
        cancelAnimationFrame(state.rafId);
        state.rafId = null;
      }
    };

    const setPaused = (value) => {
      if (state.paused === value) {
        return;
      }
      state.paused = value;
      belt.classList.toggle('logo-belt--open', value);
    };

    const toggleBeltState = () => {
      const activeBox = belt.querySelector('.materialboxed.active');
      setPaused(Boolean(activeBox));
    };

    const observer = new MutationObserver(toggleBeltState);
    track.querySelectorAll('.materialboxed').forEach((img) => {
      observer.observe(img, { attributes: true, attributeFilter: ['class'] });
    });

    const handleOverlayClick = (event) => {
      if (event.target.classList.contains('materialbox-overlay')) {
        setTimeout(toggleBeltState, 150);
      }
    };

    const handleKeyUp = (event) => {
      if (event.key === 'Escape' || event.key === 'Esc') {
        setTimeout(toggleBeltState, 150);
      }
    };

    document.addEventListener('click', handleOverlayClick);
    document.addEventListener('keyup', handleKeyUp);

    const refreshMetrics = () => {
      updateMetrics();
      state.lastTs = null;
    };

    let resizeId = null;
    window.addEventListener('resize', () => {
      if (resizeId) {
        cancelAnimationFrame(resizeId);
      }
      resizeId = requestAnimationFrame(() => {
        refreshMetrics();
        resizeId = null;
      });
    });

    const handleMotionPreferenceChange = (event) => {
      if (event.matches) {
        setPaused(true);
        stopLoop();
        track.style.transform = 'none';
      } else {
        state.offset = 0;
        applyTransform();
        refreshMetrics();
        setPaused(false);
        startLoop();
      }
    };

    if (typeof prefersReduced.addEventListener === 'function') {
      prefersReduced.addEventListener('change', handleMotionPreferenceChange);
    } else if (typeof prefersReduced.addListener === 'function') {
      prefersReduced.addListener(handleMotionPreferenceChange);
    }

    const beltObserver = new MutationObserver((mutations) => {
      if (mutations.some((mutation) => mutation.attributeName === 'style')) {
        refreshMetrics();
      }
    });
    beltObserver.observe(belt, { attributes: true, attributeFilter: ['style'] });

    refreshMetrics();
    applyTransform();
    startLoop();
  });
})();
