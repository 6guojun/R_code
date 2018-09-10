<template>
  <div :class="['vuetable-pagination-info', css.infoClass]"
       v-html="paginationInfo" class="large text">
  </div>
</template>

<script>
import PaginationInfoMixin from './VuetablePaginationInfoMixin.vue'

export default {
  mixins: [PaginationInfoMixin],
}
</script>

<style scoped>
    .large.text {
       font-size: 1.2rem;
    }
</style>