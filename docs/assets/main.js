const search = document.querySelector('#search');
const documented = document.querySelector('#documented');
function filterEntries() {
  const query = search.value.trim().toLocaleLowerCase();
  let visible = 0;
  document.querySelectorAll('.entry').forEach(entry => {
    entry.hidden = !entry.dataset.search.includes(query) ||
      (documented.checked && entry.dataset.documented !== '1');
    if (!entry.hidden) visible++;
  });
  document.querySelectorAll('.group').forEach(group => {
    group.hidden = ![...group.querySelectorAll('.entry')].some(entry => !entry.hidden);
  });
  document.querySelector('#empty').hidden = visible > 0;
}
if (search) {
  search.addEventListener('input', filterEntries);
  documented.addEventListener('change', filterEntries);
}
document.querySelectorAll('[data-copy]').forEach(button => {
  button.addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(document.getElementById(button.dataset.copy).textContent);
      button.textContent = 'コピーしました';
    } catch {
      button.textContent = 'コードを選択してコピーしてください';
    }
  });
});
